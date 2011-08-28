// Begin CVS Header
//   $Source: /fs/turing/cvs/copasi_dev/copasi/tssanalysis/CILDMMethod.cpp,v $
//   $Revision: 1.32 $
//   $Name: Build-33 $
//   $Author: shoops $
//   $Date: 2010/07/16 19:03:29 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

#include "copasi.h"
#include "CILDMMethod.h"
#include "CTSSATask.h"
#include "CTSSAProblem.h"

#include "model/CReaction.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "model/CModel.h"
#include "model/CMetab.h"
//#include "model/CState.h"
//#include "utilities/CMatrix.h"
//#include "utilities/CAnnotatedMatrix.h"
//#include "report/CCopasiObjectReference.h"

#include "clapackwrap.h"        // CLAPACK
#include "blaswrap.h"           // BLAS

//#define ILDMDEBUG

CILDMMethod::CILDMMethod(const CCopasiContainer * pParent):
    CTSSAMethod(CCopasiMethod::tssILDM, pParent) //,
    // mpState(NULL),
    // mY(NULL)
{
  //assert((void *) &mData == (void *) &mData.dim);

  // addObjectReference("Number of slow variables", mSlow, CCopasiObject::ValueInt);
  // addMatrixReference("Contribution of Metabolites to Slow Space", mVslow, CCopasiObject::ValueDbl);

  mData.pMethod = this;
  initializeParameter();
}

CILDMMethod::CILDMMethod(const CILDMMethod & src,
                         const CCopasiContainer * pParent):
    CTSSAMethod(src, pParent) //,
    //mpState(NULL),
    //mY(NULL)
{
  //assert((void *) &mData == (void *) &mData.dim);

  mData.pMethod = this;
  initializeParameter();
}

CILDMMethod::~CILDMMethod()
{
  pdelete(mpState);
}

void CILDMMethod::initializeParameter()
{
  addObjectReference("Number of slow variables", mSlow, CCopasiObject::ValueInt);
  addMatrixReference("Contribution of Species to Slow Space", mVslow, CCopasiObject::ValueDbl);

  initializeIntegrationsParameter();

  assertParameter("Deuflhard Tolerance", CCopasiParameter::UDOUBLE, (C_FLOAT64) 1.0e-6);


  createAnnotationsM();
  emptyVectors();
}

void CILDMMethod::step(const double & deltaT)
{
  C_INT failed_while = 0;

  C_INT dim = mData.dim;
  C_INT fast = 0;
  C_INT slow = dim - fast;

  C_INT i, j, k, info_schur = 0;
  mY_initial.resize(dim);
  mJacobian_initial.resize(dim, dim);
  mQ.resize(dim, dim);
  mR.resize(dim, dim);

  mQ_desc.resize(dim, dim);
  mR_desc.resize(dim, dim);

  mTd.resize(dim, dim);
  mTdInverse.resize(dim, dim);
  mQz.resize(dim, dim);

  mTd_save.resize(dim, dim);
  mTdInverse_save.resize(dim, dim);

  mpModel->updateSimulatedValues(mReducedModel);
  // TO REMOVE : mpModel->applyAssignments();
  mpModel->calculateJacobianX(mJacobian, 1e-6, 1e-12);

  C_INT flag_jacob;
  flag_jacob = 1;  // Set flag_jacob=0 to printing Jacobian

  // To get the reduced Stoichiometry Matrix;

  CMatrix<C_FLOAT64> Stoichiom;
  Stoichiom = mpModel -> getRedStoi();

  // const CCopasiVector< CReaction > & reacs = copasiModel->getReactions();
  C_INT32 reacs_size = mpModel -> getRedStoi().size();
  reacs_size = reacs_size / dim; //TODO what is this?

  /* the vector mY is the current state of the system*/

  C_FLOAT64 number2conc = 1.;  //= mpModel->getNumber2QuantityFactor()
  // / mpModel->getCompartments()[0]->getInitialValue();

  //this is an ugly hack that only makes sense if all metabs are in the same compartment
  //at the moment is is the only case the algorithm deals with

  CVector<C_FLOAT64> Xconc; //current state converted to concentrations
  Xconc.resize(dim);

  for (i = 0; i < dim; ++i)
    Xconc[i] = mY[i] * number2conc;

  for (i = 0; i < dim; i++)
    mY_initial[i] = mY[i];

  CVector<C_FLOAT64> Xconc_initial; //current state converted to concentrations
  Xconc_initial.resize(dim);

  for (i = 0; i < dim; ++i)
    Xconc_initial[i] = mY_initial[i] * number2conc;

  // save initial  Jacobian before next time step
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      mJacobian_initial(i, j) = mJacobian(i, j);

  // Next time step
  integrationStep(deltaT);

  mpModel->updateSimulatedValues(mReducedModel);
  // TO REMOVE : mpModel->applyAssignments();

  // Calculate Jacobian for time step control
  mpModel->calculateJacobianX(mJacobian, 1e-6, 1e-12);

#ifdef ILDMDEBUG

  if (flag_jacob == 0)
    {
      std::cout << "Jacobian_next:" << std::endl;
      std::cout << mJacobian << std::endl;
    }

#endif

// to be removed
  CMatrix<C_FLOAT64> Jacobian_for_test;
  Jacobian_for_test.resize(dim, dim);

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      Jacobian_for_test(i, j) = mJacobian_initial(i, j);

  // save initial  Jacobian before next time step
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      mJacobian_initial(i, j) = mJacobian(i, j);

  schur(info_schur);

#ifdef ILDMDEBUG
  std::cout << "Eigenvalues of next Jacobian" << std::endl;

  for (i = 0; i < dim; i++)
    std::cout << mR(i, i) << std::endl;

#endif

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      mJacobian_initial(i, j) = Jacobian_for_test(i, j);

  //end to be removed

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      {
        mTd_save(i, j) = 0;
        mTdInverse_save(i, j) = 0;
      }

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      {
        mTd(i, j) = 0;
        mTdInverse(i, j) = 0;
      }

  /** Schur  Decomposition of Jacobian (reordered).
  Output:  mQ - transformation matrix mR - block upper triangular matrix (with ordered eigenvalues) */

  C_INT failed = 0;
//  C_INT info_schur = 0;

  schur(info_schur); // TO DO : move the test to the TSSAMethod

  if (info_schur)
    {
      CCopasiMessage(CCopasiMessage::WARNING,
                     MCTSSAMethod + 5, mTime - deltaT);

      goto integration;
    }

  C_INT flag_schur;
  flag_schur = 0;

#ifdef ILDMDEBUG

  std::cout << "Eigenvalues of initial Jacobian" << std::endl;

  for (i = 0; i < dim; i++)
    std::cout << mR(i, i) << std::endl;

  if (flag_schur == 1)
    {
      std::cout << "Schur Decomposition" << std::endl;
      std::cout << "mR - block upper triangular matrix :" << std::endl;
      std::cout << mR << std::endl;
    }

#endif

  /* If complex eigenvalues */

  //BUG 873
  if (mR(dim - 1, dim - 1) == mR(dim - 2 , dim - 2))
    if (dim == 2)
      {
        slow = dim;
        goto integration;
      }

  // No reduction if the smallest eigenvalue is positive

  if (mR(dim - 1, dim - 1) >= 0)
    {
      slow = dim;
      fast = 0;

      CCopasiMessage(CCopasiMessage::WARNING,
                     MCTSSAMethod + 6, mTime - deltaT);

      failed = 1;
      goto integration;
    }

  // C_INT number, k;
  /** Classical ILDM iterations. The number of slow variables is decreased until the Deuflhard criterium holds */
  /*  do START slow iterations */

  while (slow > 1)
    {

      fast = fast + 1;
      slow = dim - fast;

      if (fast < dim - 1)
        if (mR(slow, slow) == mR(slow - 1 , slow - 1))
          fast = fast + 1;

      slow = dim - fast;

      for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
          {
            mTd_save(i, j) = mTd(i, j);
            mTdInverse_save(i, j) = mTdInverse(i, j);
          }

      C_INT info = 0;
      failed_while = 0;

      if (slow == 0)
        {
          for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
              {
                mTdInverse(i, j) = mQ(j, i);
                mTd(i, j) = mQ(i, j);
                mQz(i, j) = mR(i, j);
              }
        }
      else
        {
          /** Solution of Sylvester equation for given slow, mQ,mR
                Output: mTd, mTdinverse and mQz (mQz is used later for newton iterations) */
          sylvester(slow, info);

          if (info)
            {
              CCopasiMessage(CCopasiMessage::WARNING,
                             MCTSSAMethod + 7, slow, mTime - deltaT);

              failed_while = 1;
              goto integration;
            }
        }

      /* Check real parts of eigenvalues of Jacobian */

      for (i = slow ; i < dim; i++)
        if (mR(i , i) >= 0)
          {
            failed_while = 1;
            goto integration;
          }

      if (fast > 0)
        mEPS = 1 / fabs(mR(slow , slow));

      mCfast.resize(fast);

      /** Deuflhard Iteration:  Prove Deuflhard criteria, find consistent initial value for DAE
      output:  info - if Deuflhard is satisfied for this slow;
      transformation matrices mTd and mTdinverse */

      info = 0;

      C_INT help;
      help = 0;

      deuflhard(slow, info);
      help = help + 1;

      failed_while = 0;

      if (info)
        {
          failed_while = 1;
          goto integration;
        }
    }

  /** end of iterations to find the number of slow modes */

integration:

  if ((failed == 1) || (failed_while == 1))
    {
      if (slow < dim)
        {
          fast = fast - 1;
          slow = dim - fast;

          if ((fast >= 1) && (mR(slow - 1, slow - 1) == mR(slow , slow)))
            fast = fast - 1;

          slow = dim - fast;

          for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
              {
                mTd(i, j) = mTd_save(i, j);
                mTdInverse(i, j) = mTdInverse_save(i, j);
              }
        }
    }

  mSlow = slow;

  if (slow == dim)
    CCopasiMessage(CCopasiMessage::WARNING,
                   MCTSSAMethod + 8, mTime);

  //  test for orthogonality of the slow space

  C_INT flag_orthog = 1;

  if (flag_orthog == 0)
    {
      CMatrix<C_FLOAT64> orthog_prove;
      orthog_prove.resize(dim, dim);

      for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
          orthog_prove(i, j) = orthog(i, j);
    }

  C_INT flag_develop;
  flag_develop = 0; //1;

  if (flag_develop == 0)
    {
      C_INT info_schur_desc = 0;

      /** Schur  Decomposition of Jacobian (with another  sorting: from fast to slow).
        Output:  mQ_desc - transformation matrix
               mR_desc - block upper triangular matrix (with ordered eigenvalues) */

      schur_desc(info_schur_desc);

      for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
          {
            mTd_save(i, j) = mTd(i, j);
            mTdInverse_save(i, j) = mTdInverse(i, j);
          }

      for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
          {
            mTdInverse(i, j) = mQ_desc(j, i);
            mTd(i, j) = mQ_desc(i, j);
          }

      C_INT slow_desc;
      slow_desc = fast;

      mat_anal_fast_space(slow_desc);
      mat_anal_mod_space(slow_desc);

      CVector<C_FLOAT64> Reac_slow_space_orth;
      Reac_slow_space_orth.resize(reacs_size);

      for (i = 0; i < reacs_size; i ++)
        {
          Reac_slow_space_orth[i] = 0;

          for (j = 0; j < dim; j++)
            {
              Reac_slow_space_orth[i] = Reac_slow_space_orth[i] + mVfast_space[j] * Stoichiom(j, i);
            }
        }

      for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
          {
            mTd(i, j) = mTd_save(i, j);
            mTdInverse(i, j) = mTdInverse_save(i, j);
          }
    }

  // end of test

  // Post Analysis

  mat_anal_mod(slow);
  mat_anal_metab(slow);
  mat_anal_mod_space(slow);
  mat_anal_fast_space(slow);

  /** This part corresponds to the investigation of fast and slow reactions.  In  development at
   the moment. To activate the calculations take flag_develop = 0;
  */

  if (flag_develop == 0)
    {
      CMatrix<C_FLOAT64> Slow_react_contr;
      Slow_react_contr.resize(dim, reacs_size);

      CMatrix<C_FLOAT64> Mode_react_contr;
      Mode_react_contr.resize(dim, reacs_size);

      CVector<C_FLOAT64> Reac_slow_space;
      Reac_slow_space.resize(reacs_size);

      CVector<C_FLOAT64> Reac_fast_space;
      Reac_fast_space.resize(reacs_size);

      mReacSlowSpace.resize(reacs_size); //NEW TAB

      for (i = 0; i < dim; i ++)
        for (j = 0; j < reacs_size; j++)
          {
            Slow_react_contr(i, j) = 0;

            for (k = 0; k < dim; k++)
              Slow_react_contr(i, j) = Slow_react_contr(i, j) + mVslow(i, k) * Stoichiom(k, j);
          }

      CVector<C_FLOAT64> denom_mode;
      denom_mode.resize(dim);

      for (j = 0; j < dim; j++)
        denom_mode[j] = 0;

      for (i = 0; i < dim; i++)
        for (j = 0; j < reacs_size; j++)
          denom_mode[i] = denom_mode[i] + fabs(Slow_react_contr(i, j));

      for (i = 0; i < dim; i++)
        for (j = 0; j < reacs_size; j++)
          {
            if (denom_mode[i] == 0)
              Mode_react_contr(i, j) = 0;
            else
              Mode_react_contr(i, j) = (Slow_react_contr(i, j)) / denom_mode[i] * 100;
          }

      CVector<C_FLOAT64> denom_reac;
      denom_reac.resize(reacs_size);

      for (j = 0; j < reacs_size; j++)
        denom_reac[j] = 0;

      for (i = 0; i < reacs_size; i++)
        for (j = 0; j < dim; j++)
          denom_reac[i] = denom_reac[i] + fabs(Slow_react_contr(j, i));

      for (i = 0; i < reacs_size; i++)
        for (j = 0; j < dim; j++)
          {
            if (denom_reac[i] == 0)
              Slow_react_contr(j, i) = 0;
            else
              Slow_react_contr(j, i) = (Slow_react_contr(j , i)) / denom_reac[i] * 100;
          }

      for (i = 0; i < reacs_size; i ++)
        {
          Reac_slow_space[i] = 0;

          for (j = 0; j < dim; j++)
            {
              Reac_slow_space[i] = Reac_slow_space[i] + mVslow_space[j] * Stoichiom(j, i);
            }
        }

      C_FLOAT64 length;
      length = 0;

      for (i = 0; i < reacs_size; i++)
        length = length + fabs(Reac_slow_space[i]);

      for (i = 0; i < reacs_size; i++)
        if (length > 0)
          mReacSlowSpace[i] = Reac_slow_space[i] / length * 100;
        else
          mReacSlowSpace[i] = 0;

      for (i = 0; i < reacs_size; i ++)
        {
          Reac_fast_space[i] = 0;

          for (j = 0; j < dim; j++)
            Reac_fast_space[i] = Reac_fast_space[i] + mVfast_space[j] * Stoichiom(j, i);
        }

      length = 0;

      for (i = 0; i < reacs_size; i++)
        length = length + fabs(Reac_fast_space[i]);

      for (i = 0; i < reacs_size; i++)
        if (length > 0)
          Reac_fast_space[i] = Reac_fast_space[i] / length * 100;
        else
          Reac_fast_space[i] = 0;

#ifdef ILDMDEBUG

      std::cout << "**********************************************************" << std::endl;
      std::cout << "**********************************************************" << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
#endif

      mTMP1.resize(reacs_size, dim);
      mTMP2.resize(reacs_size, dim);
      mTMP3.resize(reacs_size, 1);

      for (i = 0; i < dim; i++)
        for (j = 0; j < reacs_size; j++)
          {
            mTMP1(j, i) = Mode_react_contr(i, j);
            mTMP2(j, i) = Slow_react_contr(i, j);
          }

      for (j = 0; j < reacs_size; j++)
        mTMP3(j, 0) = Reac_fast_space[j];
    }

  // End of reaction analysis

  mpModel->updateSimulatedValues(mReducedModel);
  // TO REMOVE : mpModel->applyAssignments();

  // Calculate Jacobian for time step control

  mpModel->calculateJacobianX(mJacobian, 1e-6, 1e-12);

  // new entry for every entry contains the current data of currently step
  setVectors(slow);

  // set the stepcounter
  mCurrentStep += 1;

  return;
}

/** Newton: Looking for consistent initial value for DAE system
Output:  mCfast, info
 */

void CILDMMethod::newton(C_FLOAT64 *ys, C_INT & slow, C_INT & info)
{
  C_INT i, j, iter, iterations, itermax;
  C_INT nrhs, ok, fast;

  C_FLOAT64 tol, err;
  C_INT dim = mData.dim;

  fast = dim - slow;

  CVector<C_INT> ipiv;
  ipiv.resize(fast);

  CVector<C_FLOAT64> s_22_array;
  s_22_array.resize(fast*fast);

  CVector<C_FLOAT64> gf_newton;
  gf_newton.resize(fast);

  CVector<C_FLOAT64> d_yf;
  d_yf.resize(fast);

  CVector<C_FLOAT64> y_newton;
  y_newton.resize(dim);

  CVector<C_FLOAT64> yf_newton;
  yf_newton.resize(fast);

  CVector<C_FLOAT64> x_newton;
  x_newton.resize(dim);

  CVector<C_FLOAT64> dxdt_newton;
  dxdt_newton.resize(dim);

  CVector<C_FLOAT64> g_newton;
  g_newton.resize(dim);

  CMatrix<C_FLOAT64> S_22;
  S_22.resize(fast, fast);

  C_FLOAT64 g1, g2 = 0;

  nrhs = 1;
  //tol = 1e-9;
  tol = 1e-9 / mpModel->getNumber2QuantityFactor();


  err = 10.0 / mpModel->getNumber2QuantityFactor();
  iter = 0;

  itermax = 100;
  iterations = 0;

  info = 0;

  for (i = 0; i < fast; i++)
    for (j = 0; j < fast; j++)
      S_22(i, j) = mQz(i, j);

  for (i = 0; i < fast; i++)
    yf_newton[i] = mCfast[i];

  for (i = 0; i < fast; i++)
    for (j = 0; j < fast; j++)
      s_22_array[j + fast*i] = S_22(j, i);

  for (i = 0; i < fast; i++)
    d_yf[i] = 0.;

  while (err > tol)
    {
      iter ++;

      if (iter > itermax)
        {

          info = 1;
          return;
        }

      for (i = 0; i < fast; i++)
        yf_newton[i] = yf_newton[i] + d_yf[i];

      /* back transformation */

      for (i = 0; i < slow; i++)
        y_newton[i] = ys[i];

      for (i = slow; i < dim; i++)
        y_newton[i] = yf_newton[i - slow];

      for (i = 0; i < dim; i++)
        {
          x_newton[i] = 0.0;

          for (j = 0; j < dim; j++)
            x_newton[i] = x_newton[i] + mTd(i, j) * y_newton[j];
        }

      calculateDerivativesX(x_newton.array(), dxdt_newton.array());

      for (i = 0; i < dim; i++)
        {
          g_newton[i] = 0.;

          for (j = 0; j < dim; j++)
            g_newton[i] = g_newton[i] + mTdInverse(i, j) * dxdt_newton[j];
        }

      // for (i = 0; i < fast; i++)
      //  gf_newton[i] = -1. * g_newton[i + slow];

      for (i = 0; i < fast; i++)
        {
          gf_newton[i] = -1. * g_newton[i + slow];

        }

      /*       int dgesv_(integer *n, integer *nrhs, doublereal *a, integer
       * *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info)
       *
       *  -- LAPACK driver routine (version 3.0) --
       *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
       *     Courant Institute, Argonne National Lab, and Rice University
       *     March 31, 1993
       *
       *
       *  Purpose
       *  =======
       *
       *  DGESV computes the solution to a real system of linear equations
       *     A * X = B,
       *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
       *
       *  The LU decomposition with partial pivoting and row interchanges is
       *  used to factor A as
       *     A = P * L * U,
       *  where P is a permutation matrix, L is unit lower triangular, and U is
       *  upper triangular.  The factored form of A is then used to solve the
       *  system of equations A * X = B.
       *
       *  Arguments
       *  =========
       *
       *  N       (input) INTEGER
       *          The number of linear equations, i.e., the order of the
       *          matrix A.  N >= 0.
       *
       *  NRHS    (input) INTEGER
       *          The number of right hand sides, i.e., the number of columns
       *          of the matrix B.  NRHS >= 0.
       *
       *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
       *          On entry, the N-by-N coefficient matrix A.
       *          On exit, the factors L and U from the factorization
       *          A = P*L*U; the unit diagonal elements of L are not stored.
       *
       *  LDA     (input) INTEGER
       *          The leading dimension of the array A.  LDA >= max(1,N).
       *
       *  IPIV    (output) INTEGER array, dimension (N)
       *          The pivot indices that define the permutation matrix P;
       *          row i of the matrix was interchanged with row IPIV(i).
       *
       *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
       *          On entry, the N-by-NRHS matrix of right hand side matrix B.
       *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
       *
       *  LDB     (input) INTEGER
       *          The leading dimension of the array B.  LDB >= max(1,N).
       *
       *
       * INFO    (output) INTEGER
       *          = 0:  successful exit
       *          < 0:  if INFO = -i, the i-th argument had an illegal value
       *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
       *                has been completed, but the factor U is exactly
       *                singular, so the solution could not be computed.
       */

      dgesv_(&fast, &nrhs, s_22_array.array(), &fast, ipiv.array(), gf_newton.array(), &fast, &ok);

      if (ok != 0)
        {

          info = 2;
          break;
        }

      for (i = 0; i < fast; i++)
        d_yf[i] = gf_newton[i];

      err = -10.;

      /*      for (i = 0; i < fast; i++)
              {
                gf_newton[i] = fabs(gf_newton[i]);

                if (err < gf_newton[i])
                  err = gf_newton[i];
              }

      */


      for (i = 0; i < fast; i++)
        {
          gf_newton[i] = fabs(gf_newton[i]);

          if (err < gf_newton[i])
            err = gf_newton[i];
        }




      iterations = iterations + 1;

      /* stop criterion of newton method */

      //  C_FLOAT64 g1, g2;

      //  g2 = err;

      if (iter == 1)
        g1 = 3.0 * err;
      else
        g1 = g2;

      g2 = err;


      if (g2 / g1 > 1.0)
        {

          info = 1;
          break;
        }
    } /* end while */

  for (i = 0; i < fast; i++)
    mCfast[i] = yf_newton[i];

  return;
}

void CILDMMethod::start(const CState * initialState)
{

  integrationMethodStart(initialState);

  mDtol = * getValue("Deuflhard Tolerance").pUDOUBLE;

  /* ILDM related staff  */

  //mDtol = mpProblem->getDeufelhardTol();

  mVslow.resize(mData.dim, mData.dim);
  mVslow_metab.resize(mData.dim, mData.dim);
  mVslow_space.resize(mData.dim);
  mVfast_space.resize(mData.dim);

  //createAnnotationsM();
  emptyVectors();

  return;
}

/**
  Deuflhard Iteration:  Prove Deuflhard criteria, find consistent initial value for DAE
  output:  info - if Deuflhard is satisfied for given slow; transformation matrices
  mTd and mTdinverse
 */
void CILDMMethod::deuflhard(C_INT & slow, C_INT & info)
{
  C_INT i, j;
  C_INT dim = mData.dim;
  C_INT fast = dim - slow;
  C_INT flag_deufl;


  flag_deufl = 1;  // set flag_deufl = 0 to print the results of calculations

  /* calculations before relaxing yf to slow manifold */

  CVector<C_FLOAT64> c_full;
  c_full.resize(dim);

  CVector<C_FLOAT64> c_slow;
  c_slow.resize(slow);

  /* the vector mY is the current state of the system*/

  C_FLOAT64 number2conc =  1.; // mpModel->getNumber2QuantityFactor()
  // / mpModel->getCompartments()[0]->getInitialValue();

  //this is an ugly hack that only makes sense if all metabs are in the same compartment
  //at the moment is is the only case the algorithm deals with

  CVector<C_FLOAT64> Xconc; //current state converted to concentrations
  Xconc.resize(dim);

  for (i = 0; i < dim; ++i)
    Xconc[i] = mY_initial[i] * number2conc;

  for (i = 0; i < dim; i++)
    {
      c_full[i] = 0.0;

      for (j = 0; j < dim; j++)
        c_full[i] = c_full[i] + mTdInverse(i, j) * Xconc[j];
    }

  for (j = 0; j < slow; j++)
    c_slow[j] = c_full[j];

  for (j = 0; j < fast; j++)
    mCfast[j] = c_full[j + slow];

  CVector<C_FLOAT64> g_full;
  g_full.resize(dim);

  CVector<C_FLOAT64> g_slow;
  g_slow.resize(slow);

  CVector<C_FLOAT64> g_fast;
  g_fast.resize(fast);

  CVector<C_FLOAT64> dxdt;
  dxdt.resize(dim);

  mpModel->updateSimulatedValues(mReducedModel);

  for (j = 0; j < dim; j++)
    dxdt[j] = 0.;

  CVector<C_FLOAT64> x_help;
  x_help.resize(dim);

  for (j = 0; j < dim; j++)
    {
      x_help[j] = mY_initial[j] * number2conc;
    }

  calculateDerivativesX(x_help.array(), dxdt.array());

  for (i = 0; i < dim; i++)
    {
      g_full[i] = 0.0;

      for (j = 0; j < dim; j++)
        g_full[i] = g_full[i] + mTdInverse(i, j) * dxdt[j];
    }

  for (j = 0; j < slow; j++)
    g_slow[j] = g_full[j];

  info = 0;

  /**  NEWTON: Looking for consistent initial value for DAE system
  Output:  mCfast, info */
  newton(c_slow.array(), slow, info);


  if (info)
    {
      /* TODO */

      return;
    }


  /* calculation of g_relax at point x_relax (after relaxing yf to slow manifold)*/

  CVector<C_FLOAT64> c_relax;
  c_relax.resize(dim);

  CVector<C_FLOAT64> x_relax;
  x_relax.resize(dim);

  CVector<C_FLOAT64> dxdt_relax;
  dxdt_relax.resize(dim);

  CVector<C_FLOAT64> g_relax;
  g_relax.resize(dim);

  for (i = 0; i < slow; i++)
    c_relax[i] = c_slow[i];

  for (i = slow; i < dim; i++)
    c_relax[i] = mCfast[i - slow];

  for (i = 0; i < dim; i++)
    {
      x_relax[i] = 0.0;

      for (j = 0; j < dim; j++)
        x_relax[i] = x_relax[i] + mTd(i, j) * c_relax[j];
    }

  calculateDerivativesX(x_relax.array(), dxdt_relax.array());

  for (i = 0; i < dim; i++)
    {
      g_relax[i] = 0.0;

      for (j = 0; j < dim; j++)
        g_relax[i] = g_relax[i] + mTdInverse(i, j) * dxdt_relax[j];
    }

  CVector<C_FLOAT64> re;
  re.resize(slow);

  /* stop criterion for slow reaction modes */

  for (i = 0; i < slow; i++)
    {
      re[i] = fabs(g_relax[i] - g_slow[i]);
      re[i] = re[i] * mEPS;



    }

  C_FLOAT64 max = 0.;

  for (i = 0; i < slow; i++)
    if (max < re[i])
      max = re[i];

  C_FLOAT64 max1;
  C_FLOAT64 norm = 0;

  for (i = 0; i < slow; i++)
    norm = norm + fabs(g_relax[i] - g_slow[i]);

  max1 = norm * mEPS;

  if (max >= mDtol / mpModel->getNumber2QuantityFactor())
    info = 1;
  else
    info = 0;

  return;
}

/**
 * Empty every vector to be able to fill them with new values for a new calculation.
 * Also nullify the step counter.
 **/
void CILDMMethod::emptyVectors()
{
  mCurrentStep = 0;
  mVec_mVslow.erase(mVec_mVslow.begin(), mVec_mVslow.end());
  mVec_TimeScale.erase(mVec_TimeScale.begin(), mVec_TimeScale.end());
  mVec_mVslowMetab.erase(mVec_mVslowMetab.begin(), mVec_mVslowMetab.end());
  mVec_mVslowSpace.erase(mVec_mVslowSpace.begin(), mVec_mVslowSpace.end());
  mVec_SlowModes.erase(mVec_SlowModes.begin(), mVec_SlowModes.end());

  /* temporary tabs */

  mVec_mTMP1.erase(mVec_mTMP1.begin(), mVec_mTMP1.end());
  mVec_mTMP2.erase(mVec_mTMP2.begin(), mVec_mTMP2.end());
  mVec_mTMP3.erase(mVec_mTMP3.begin(), mVec_mTMP3.end());


}

/**
 *upgrade all vectors with values from actually calculalion for current step
 **/
void CILDMMethod::setVectors(int slowMode)
{
  mVec_mVslow.push_back(mCurrentStep);
  mVec_mVslow[mCurrentStep].resize(mData.dim, mData.dim);
  mVec_mVslow[mCurrentStep] = mVslow;

  mVec_TimeScale.push_back(mCurrentStep);
  mVec_TimeScale[mCurrentStep].resize(mData.dim);
  unsigned C_INT32 i;

  for (i = 0; i < (unsigned C_INT32) mData.dim; i++)
    mVec_TimeScale[mCurrentStep][i] = -1 / mR(i, i);

  mVec_mVslowMetab.push_back(mCurrentStep);
  mVec_mVslowMetab[mCurrentStep].resize(mData.dim, mData.dim);
  mVec_mVslowMetab[mCurrentStep] = mVslow_metab;

  mVec_mVslowSpace.push_back(mCurrentStep);
  mVec_mVslowSpace[mCurrentStep].resize(mData.dim);

  for (i = 0; i < mVslow_space.size(); i++)
    {
      mVec_mVslowSpace[mCurrentStep][i] = mVslow_space[i];
    }

  mVec_mVfastSpace.push_back(mCurrentStep);
  mVec_mVfastSpace[mCurrentStep].resize(mData.dim);
  mVec_mVfastSpace[mCurrentStep] = mVfast_space;

  mVec_SlowModes.push_back(mCurrentStep);
  mVec_SlowModes[mCurrentStep] = slowMode;

  mCurrentTime.push_back(mCurrentStep);
  mCurrentTime[mCurrentStep] = mTime;

  // NEW TAB

  mVec_mReacSlowSpace.push_back(mCurrentStep);
  mVec_mReacSlowSpace[mCurrentStep].resize(mReacSlowSpace.size());
  mVec_mReacSlowSpace[mCurrentStep] = mReacSlowSpace;

  /* temporary tabs */

  C_INT reacs_size = mpModel->getReactions().size();

  mVec_mTMP1.push_back(mCurrentStep);
  mVec_mTMP1[mCurrentStep].resize(reacs_size, mData.dim);
  mVec_mTMP1[mCurrentStep] = mTMP1;

  mVec_mTMP2.push_back(mCurrentStep);
  mVec_mTMP2[mCurrentStep].resize(reacs_size, mData.dim);
  mVec_mTMP2[mCurrentStep] = mTMP2;

  mVec_mTMP3.push_back(mCurrentStep);
  mVec_mTMP3[mCurrentStep].resize(reacs_size, 1);
  mVec_mTMP3[mCurrentStep] = mTMP3;
}
/**
 * Create the CArraAnnotations for every ILDM-tab in the CQTSSAResultSubWidget.
 * Input for each CArraAnnotations is a seperate CMatrix.
 **/
void CILDMMethod::createAnnotationsM()
{
  tableNames.erase(tableNames.begin(), tableNames.end());

  std::string name;

  name = "Contribution of species to modes";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTmp1 = new CArrayAnnotation("Contribution of species to modes", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mVslowPrint), true);
  pTmp1->setMode(0, pTmp1->STRINGS);
  pTmp1->setMode(1, pTmp1->VECTOR);
  pTmp1->setDescription(" ");
  pTmp1->setDimensionDescription(0, "Contribution to  mode (TS - corresponding timescale)");
  pTmp1->setDimensionDescription(1, "Species");
  pVslowPrintAnn = pTmp1;

  mapTableToName[name] = pVslowPrintAnn;

  name = "Modes distribution for species";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTmp2 = new CArrayAnnotation("Modes distribution for species", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mVslowMetabPrint), true);
  pTmp2->setMode(1, pTmp2->STRINGS);
  pTmp2->setMode(0, pTmp2->VECTOR);
  pTmp2->setDescription(" ");
  pTmp2->setDimensionDescription(0, "Mode distribution for each species");
  pTmp2->setDimensionDescription(1, "Modes (TS - corresponding  timescale)");
  pVslowMetabPrintAnn = pTmp2;

  mapTableToName[name] = pVslowMetabPrintAnn;

  name = "Slow space";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTmp3 = new CArrayAnnotation("Slow space", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mVslowSpacePrint), true);
  pTmp3->setMode(1, pTmp3->STRINGS);
  pTmp3->setMode(0, pTmp3->VECTOR);
  pTmp3->setDescription(" ");
  pTmp3->setDimensionDescription(0, "Species");
  pTmp3->setDimensionDescription(1, "Contribution to slow space");
  pVslowSpacePrintAnn = pTmp3;

  mapTableToName[name] = pVslowSpacePrintAnn;

  name = "Fast space";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTmp4 = new CArrayAnnotation("Fast space", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mVfastSpacePrint), true);
  pTmp4->setMode(1, pTmp4->STRINGS);
  pTmp4->setMode(0, pTmp4->VECTOR);
  pTmp4->setDescription(" ");
  pTmp4->setDimensionDescription(0, "Species");
  pTmp4->setDimensionDescription(1, "Contribution to fast space");
  pVfastSpacePrintAnn = pTmp4;

  mapTableToName[name] = pVfastSpacePrintAnn;

  // NEW TAB

  name = "Reactions slow space";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTmp5 = new CArrayAnnotation("Reactions slow space", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mReacSlowSpacePrint), true);
  pTmp5->setMode(1, pTmp5->STRINGS);
  pTmp5->setMode(0, pTmp5->VECTOR);
  pTmp5->setDescription(" ");
  pTmp5->setDimensionDescription(0, "Reactions");
  pTmp5->setDimensionDescription(1, "Contribution to slow space");
  pReacSlowSpacePrintAnn = pTmp5;

  mapTableToName[name] = pReacSlowSpacePrintAnn;

  /* tamporary tabs */

  name = "Reactions contribution to the mode";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTMP1 = new CArrayAnnotation("", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mTMP1Print), true);
  pTMP1->setMode(0, pTMP1->VECTOR);
  pTMP1->setMode(1, pTMP1->STRINGS);
  pTMP1->setDescription("Reactions contribution to the mode ");
  pTMP1->setDimensionDescription(0, "Reactions");
  pTMP1->setDimensionDescription(1, "Modes (TS - corresponding  timescale)");
  pTMP1PrintAnn = pTMP1;

  mapTableToName[name] = pTMP1PrintAnn;

  name = "Reactions distribution between modes ";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTMP2 = new CArrayAnnotation("", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mTMP2Print), true);
  pTMP2->setMode(0, pTMP2->VECTOR);
  pTMP2->setMode(1, pTMP2->STRINGS);
  pTMP2->setDescription("Reactions distribution between modes ");
  pTMP2->setDimensionDescription(0, "Reactions");
  pTMP2->setDimensionDescription(1, "Modes (TS - corresponding  timescale)");
  pTMP2PrintAnn = pTMP2;

  mapTableToName[name] = pTMP2PrintAnn;

  name = "Reactions fast space";
  tableNames.push_back(name);

  CArrayAnnotation *
  pTMP3 = new CArrayAnnotation("", this,
                               new CCopasiMatrixInterface<CMatrix<C_FLOAT64> >(&mTMP3Print), true);
  pTMP3->setMode(0, pTMP3->VECTOR);
  pTMP3->setMode(1, pTMP3->STRINGS);
  pTMP3->setDescription("Reactions fast space ");
  pTMP3->setDimensionDescription(0, "Reactions");
  pTMP3->setDimensionDescription(1, "");
  pTMP3PrintAnn = pTMP3;

  mapTableToName[name] = pTMP3PrintAnn;
}
/**
 * Set the every CArrayAnnotation for the requested step.
 * Set also the desription of CArayAnnotation for both dimensions:
 *    - dimension description could consists of some std::srings
 *      some strings contain the Time Scale values for requested step
 *    - dimension description could consists of arrays of CommonNames
 **/
//void CILDMMethod::setAnnotationM(int step)
bool CILDMMethod::setAnnotationM(int step)
{
  if (step == 0) return false;

  if (mVec_mVslow.size() == 0) return false;

  step -= 1;
  double timeScale;
  std::string str;
  std::stringstream sstr;
  sstr.str("");
  sstr.clear();
  C_INT32 i;

  mVslowPrint.resize(mData.dim, mData.dim);
  mVslowPrint = mVec_mVslow[step];
  pVslowPrintAnn->resize();
  pVslowPrintAnn->setCopasiVector(1, &mpModel->getMetabolitesX());

  for (i = 0; i < mData.dim; i++)
    {
      timeScale = mVec_TimeScale[step][i];

      if (i < mVec_SlowModes[step])
        sstr << "Slow: ";
      else
        sstr << "Fast: ";

      sstr << timeScale;
      str = sstr.str();
      pVslowPrintAnn->setAnnotationString(0, i, str);
      sstr.str("");
      sstr.clear();
    }

  mVslowMetabPrint.resize(mData.dim, mData.dim);
  mVslowMetabPrint = mVec_mVslowMetab[step];
  pVslowMetabPrintAnn->resize();
  pVslowMetabPrintAnn->setCopasiVector(0, &mpModel->getMetabolitesX());

  for (i = 0; i < mData.dim; i++)
    {
      timeScale = mVec_TimeScale[step][i];

      if (i < mVec_SlowModes[step])
        sstr << "Slow: ";
      else
        sstr << "Fast: ";

      sstr << timeScale;
      str = sstr.str();
      pVslowMetabPrintAnn->setAnnotationString(1, i, str);
      sstr.str("");
      sstr.clear();
    }

  sstr << mVec_SlowModes[step];
  // if (mVec_SlowModes[step] > 1)
  // sstr << " slow modes";
  //else
  // sstr << " slow mode";
  sstr << " slow; ";

  C_INT dim = mData.dim;
  sstr << dim - mVec_SlowModes[step];
  sstr << " fast";

  str = sstr.str();
  mVslowSpacePrint.resize(mData.dim, 1);

  for (i = 0; i < mData.dim; i++)
    mVslowSpacePrint(i, 0) = mVec_mVslowSpace[step][i];

  pVslowSpacePrintAnn->resize();
  pVslowSpacePrintAnn->setCopasiVector(0, &mpModel->getMetabolitesX());
  pVslowSpacePrintAnn->setAnnotationString(1, 0, str);

  mVfastSpacePrint.resize(mData.dim, 1);

  for (i = 0; i < mData.dim; i++)
    mVfastSpacePrint(i, 0) = mVec_mVfastSpace[step][i];

  pVfastSpacePrintAnn->resize();
  pVfastSpacePrintAnn->setCopasiVector(0, &mpModel->getMetabolitesX());
  pVfastSpacePrintAnn->setAnnotationString(1, 0, str);

  //sstr.clear();
  sstr.str("");
  str = sstr.str();

  mReacSlowSpacePrint.resize(mReacSlowSpace.size(), 1);

  for (i = 0; i < (C_INT32) mReacSlowSpace.size(); i++)
    mReacSlowSpacePrint(i, 0) = mVec_mReacSlowSpace[step][i];

  pReacSlowSpacePrintAnn->resize();
  pReacSlowSpacePrintAnn->setCopasiVector(0, &mpModel->getReactions());
  pReacSlowSpacePrintAnn->setAnnotationString(1, 0, str);

  sstr.str("");
  str = sstr.str();

  /* temporary tabs */

  C_INT reacs_size = mpModel->getReactions().size();

  mTMP1Print.resize(reacs_size, mData.dim);
  mTMP1Print = mVec_mTMP1[step];

  pTMP1PrintAnn->resize();

  for (i = 0; i < (C_INT) mVec_mTMP1[step].numCols(); i++)
    {
      timeScale = mVec_TimeScale[step][i];

      if (i < (C_INT) mVec_SlowModes[step])
        sstr << "Slow: ";
      else
        sstr << "Fast: ";

      sstr << timeScale;
      str = sstr.str();
      pTMP1PrintAnn->setAnnotationString(1, i, str);
      sstr.str("");
      sstr.clear();
    }

  pTMP1PrintAnn->setCopasiVector(0, &mpModel->getReactions());

  mTMP2Print.resize(reacs_size, mData.dim);
  mTMP2Print = mVec_mTMP2[step];

  pTMP2PrintAnn->resize();

  for (i = 0; i < (C_INT) mVec_mTMP2[step].numCols(); i++)
    {
      timeScale = mVec_TimeScale[step][i];

      if (i < (C_INT) mVec_SlowModes[step])
        sstr << "Slow: ";
      else
        sstr << "Fast: ";

      sstr << timeScale;
      str = sstr.str();
      pTMP2PrintAnn->setAnnotationString(1, i, str);
      sstr.str("");
      sstr.clear();
    }

  pTMP2PrintAnn->setCopasiVector(0, &mpModel->getReactions());

  mTMP3Print.resize(reacs_size, 1);
  mTMP3Print = mVec_mTMP3[step];
  pTMP3PrintAnn->resize();
  pTMP3PrintAnn->setCopasiVector(0, &mpModel->getReactions());
  pTMP3PrintAnn->setAnnotationString(1, 0, str);

  return true;
}

void CILDMMethod::printResult(std::ostream * ostream) const // temporary tabs are not included in Report !!!
{
  std::ostream & os = *ostream;
  double timeScale;
  C_INT i, j, istep = 0;

  C_INT32 stepNumber;

  this->print(&os);

  assert(CCopasiRootContainer::getDatamodelList()->size() > 0);
  CTSSATask* pTask =
    dynamic_cast<CTSSATask *>((*(*CCopasiRootContainer::getDatamodelList())[0]->getTaskList())["Time Scale Separation Analysis"]);

  CTSSAProblem* pProblem = dynamic_cast<CTSSAProblem*>(pTask->getProblem());

  stepNumber = pProblem->getStepNumber();

  for (istep = 0; istep < stepNumber; istep++)
    {

      os << std::endl;
      os << "****************  Time step " << istep + 1 << " **************************  " << std::endl;

      os << std::endl;

      os << "Contribution of species to modes" << std::endl;

      os << "Rows : contribution to  mode (TS - corresponding timescale)" << std::endl;
      os << "Columns: species  ";

      for (j = 0; j < mData.dim; j++)
        {
          os << mpModel->getMetabolitesX()[j]->getObjectName() << "   ";
        }

      os << std::endl;

      for (i = 0; i < mData.dim; i++)
        {
          timeScale = mVec_TimeScale[istep][i];

          if (i < mVec_SlowModes[istep])
            os << "Slow (";
          else
            os << "Fast (";

          os << timeScale << "): ";

          for (j = 0; j < mData.dim; j++)
            os << mVec_mVslow[istep][i][j] << " ";

          os << std::endl;
        }

      os << std::endl;

      os << "Modes distribution for species" << std::endl;

      os << "Rows : Mode distribution for each species" << std::endl;
      os << "Columns: Modes (TS - corresponding  timescale) ";
      os << std::endl;

      for (i = 0; i < mData.dim; i++)
        {
          timeScale = mVec_TimeScale[istep][i];

          if (i < mVec_SlowModes[istep])
            os << "Slow (";
          else
            os << "Fast (";

          os << timeScale << ")  ";
        }

      os << std::endl;

      for (j = 0; j < mData.dim; j++)
        {
          os << mpModel->getMetabolitesX()[j]->getObjectName() << "  ";

          for (i = 0; i < mData.dim; i++)
            os << mVec_mVslowMetab[istep][j][i] << "  ";

          os << std::endl;
        }

      os << std::endl;

      os << "Slow space" << std::endl;

      os << "Rows : Species" << std::endl;
      os << "Column: Contribution to slow space ";
      os << std::endl;

      os << mVec_SlowModes[istep];
      os << " slow; ";

      os << mData.dim - mVec_SlowModes[istep];
      os << " fast";
      os << std::endl;

      for (j = 0; j < mData.dim; j++)
        {
          os << mpModel->getMetabolitesX()[j]->getObjectName() << "  ";
          os << mVec_mVslowSpace[istep][j] << "  ";

          os << std::endl;
        }

      os << std::endl;
      os << "Fast space" << std::endl;

      os << "Rows : Species" << std::endl;
      os << "Column: Contribution to fast space ";
      os << std::endl;

      os << mVec_SlowModes[istep];
      os << " slow; ";

      os << mData.dim - mVec_SlowModes[istep];
      os << " fast";
      os << std::endl;

      for (j = 0; j < mData.dim; j++)
        {
          os << mpModel->getMetabolitesX()[j]->getObjectName() << "  ";
          os << mVec_mVfastSpace[istep][j] << "  ";

          os << std::endl;
        }

      os << std::endl;
      os << "Reactions slow space" << std::endl;

      os << "Rows : Reactions" << std::endl;
      os << "Column: Contribution to slow space ";
      os << std::endl;

      for (j = 0; j < (C_INT32) mpModel->getReactions().size(); j++)
        {
          os << mpModel->getReactions()[j]->getObjectName() << "  ";
          os << mVec_mReacSlowSpace[istep][j] << "  ";

          os << std::endl;
        }

      os << std::endl;
    }

  return;
}
