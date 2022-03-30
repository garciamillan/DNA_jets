#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>


struct data {
  size_t n;
  double * t;
  double * y;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;
  double *y = ((struct data *)data)->y;

    double A = gsl_vector_get (x, 0);
    double sigmasq = gsl_vector_get (x, 1);
    double mu = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
        /* Model Yi = A * exp(-0.5 * (t_i - mu) * (t_i - mu) / sigmasq ) */
      double Yi = A * exp(-0.5 * (t[i] - mu) * (t[i] - mu) / sigmasq );
      gsl_vector_set (f, i, Yi - y[i]);
    }

  return GSL_SUCCESS;
}


int expb_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;

    double A = gsl_vector_get (x, 0);
    double sigmasq = gsl_vector_get (x, 1);
    double mu = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-0.5 * (t_i - mu) * (t_i - mu) / sigmasq )  */
      /* and the xj are the parameters (A,sigmasq,mu) */
        double e = exp(-0.5 * (t[i] - mu) * (t[i] - mu) / sigmasq );
        gsl_matrix_set (J, i, 0, e );
        gsl_matrix_set (J, i, 1, 0.5 * A * ((t[i] - mu) / sigmasq) * ((t[i] - mu) / sigmasq) * e);
        gsl_matrix_set (J, i, 2, A * ((t[i] - mu) / sigmasq) * e );
    }

  return GSL_SUCCESS;
}

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr, "iter %2zu: A = %.4f, sigmasq = %.4f, mu = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          1.0 / rcond,
          gsl_blas_dnrm2(f));
}

int rgm_gsl_Gauss_fit01 (double *t, double *y, double *weights, const size_t n)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  const size_t p = 3;

  gsl_vector *f;
  gsl_matrix *J;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  struct data d = { n, t, y };
  double x_init[3] = { 1.0, 1.0, 0.0 }; /* starting values */
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  gsl_vector_view wts = gsl_vector_view_array(weights, n);
  gsl_rng * r;
  double chisq, chisq0;
  int status, info;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  /* define the function to be minimized */
  fdf.f = expb_f;
  fdf.df = expb_df;   /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = n;
  fdf.p = p;
  fdf.params = &d;

  /* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

  /* compute initial cost function */
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  /* solve the system with a maximum of 100 iterations */
  status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, callback, NULL, &info, w);

  /* compute covariance of best fit parameters */
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar (J, 0.0, covar);

  /* compute final cost */
  gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  fprintf(stdout, "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
  fprintf(stdout, "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
  fprintf(stdout, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stdout, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stdout, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
  fprintf(stdout, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stdout, "final   |f(x)| = %f\n", sqrt(chisq));

  {
    double dof = n - p;
    double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

    fprintf(stdout, "chisq/dof = %g\n", chisq / dof);

    fprintf (stdout, "A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    fprintf (stdout, "sigmasq = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    fprintf (stdout, "mu      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
  }

  fprintf (stdout, "status = %s\n", gsl_strerror (status));

  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
  gsl_rng_free (r);
    fprintf (stdout, "c'est fini\n");

  return 0;
}

void rgm_gsl_Gauss_fit (double *t, double *y, double *weights, const size_t n, double *est_A, double *est_mu, double *est_sigmasq, double *std_A, double *std_mu, double *std_sigmasq)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  const size_t p = 3;

  gsl_vector *f;
  gsl_matrix *J;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  struct data d = { n, t, y };
  double x_init[3] = { 1.0, 1.0, 0.0 }; /* starting values */
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  gsl_vector_view wts = gsl_vector_view_array(weights, n);
  gsl_rng * r;
  double chisq, chisq0;
  int status, info;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

    fprintf(stderr, "#Info: Fitting data\n");
    int i;
    for(i=0;i<n;i++) fprintf(stderr, "%d %.5f %.5f %.5f\n",i,t[i],y[i],sqrt(1./weights[i]));
    
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  /* define the function to be minimized */
  fdf.f = expb_f;
  fdf.df = expb_df;   /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = n;
  fdf.p = p;
  fdf.params = &d;

  /* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

  /* compute initial cost function */
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  /* solve the system with a maximum of 100 iterations */
  status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, callback, NULL, &info, w);

  /* compute covariance of best fit parameters */
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar (J, 0.0, covar);

  /* compute final cost */
  gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  fprintf(stderr, "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

  
    double dof = n - p;
    double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

    fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

    fprintf (stderr, "A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    fprintf (stderr, "sigmasq = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    fprintf (stderr, "mu      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
  

  fprintf (stderr, "status = %s\n", gsl_strerror (status));

  *est_A=FIT(0);
  *est_sigmasq=FIT(1);
  *est_mu=FIT(2);
    
  *std_A=c*ERR(0);
  *std_mu=c*ERR(1);
  *std_sigmasq=c*ERR(2);
  
  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
  gsl_rng_free (r);
    
    fprintf (stderr, "c'est fini\n");
}
