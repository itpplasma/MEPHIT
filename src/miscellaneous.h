//
// Created by michi on 02.02.22.
//

#ifndef PROJECT_MATHE_MASTER_MISCELLANEOUS_H
#define PROJECT_MATHE_MASTER_MISCELLANEOUS_H

#include "mfem.hpp"
#include <vector>

namespace utility::mfem_helper{

    //assume that f : mfem::Vector -> double
    template <class F_t>
    mfem::Vector get_vector_representation_of_function(mfem::FiniteElementSpace& X, F_t & f){
        mfem::Vector v(X.GetVSize());
        //define Gridfunctions
        mfem::GridFunction f_gf;
        f_gf.MakeRef(&X, v);

        //fill v with function values
        mfem::FunctionCoefficient f_coeff_fun(f);
        f_gf.ProjectCoefficient(f_coeff_fun);

        return v;
    }

    //assume that f : mfem::Vector -> double
    template<class F_t>
    double get_l2_error_of_solution(mfem::Vector const& sol, mfem::FiniteElementSpace& U, F_t& sol_analytical){
        //get appropriate integration rule:
        const int order = U.GetMaxElementOrder();
        int order_quad = std::max(2, 2*order+1);
        const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
        for (int i=0; i < mfem::Geometry::NumGeom; ++i)
        {
            irs[i] = &(mfem::IntRules.Get(i, order_quad));
        }

        //copy vector just to be sure
        mfem::Vector sol_t(sol);
        mfem::GridFunction sol_gf;
        sol_gf.MakeRef(&U, sol_t);

        //compute L2 error
        mfem::FunctionCoefficient sol_analytical_fc(sol_analytical);
        double err_L2  = sol_gf.ComputeL2Error(sol_analytical_fc, irs);

        return err_L2;
    }

    //assume that f : mfem::Vector -> double, and grad_analytical void(mfem::Vector x, mfem::Vector grad)
    template<class F_t, class F_grad_t>
    double get_h1_error_of_solution(mfem::Vector const& sol, mfem::FiniteElementSpace& U, F_t& sol_analytical, F_grad_t& grad_analytical){
        //get appropriate integration rule:
        const int order = U.GetMaxElementOrder();
        int order_quad = std::max(2, 2*order+1);
        const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
        for (int i=0; i < mfem::Geometry::NumGeom; ++i)
        {
            irs[i] = &(mfem::IntRules.Get(i, order_quad));
        }

        //copy vector just to be sure
        mfem::Vector sol_t(sol);
        mfem::GridFunction sol_gf;
        sol_gf.MakeRef(&U, sol_t);

        //compute L2 error
        mfem::FunctionCoefficient sol_analytical_fc(sol_analytical);
        const int dim = U.GetMesh()->Dimension();
        mfem::VectorFunctionCoefficient grad_analytical_fc(dim, grad_analytical);

        double const  err_H1  = sol_gf.ComputeH1Error(&sol_analytical_fc, &grad_analytical_fc, irs);

        return err_H1;
    }

    //evaluate (Au,v)
    double evaluate_bilinear_form(mfem::SparseMatrix const& A, mfem::Vector const& u, mfem::Vector const& v);
    //evaluate (Au,v)
    double evaluate_bilinear_form(mfem::Operator*  A, mfem::Vector const& u, mfem::Vector const& v);

    mfem::Vector get_random_vector(const int n);

    double get_operator_application_difference(mfem::Operator* A, mfem::Operator* B);

    mfem::SparseMatrix block_to_sparse_matrix(mfem::BlockMatrix const& B);

    mfem::Vector solve_penalty(mfem::SparseMatrix const& A, mfem::Vector const& rhs, mfem::Array<int> const& essential_dofs, mfem::Vector const& u, const double penalty_fac);


}

namespace utility::matlab_like_features{
    template <typename T>
    std::vector<T> linspace(T a, T b, size_t N) {
        T h = (b - a) / static_cast<T>(N-1);
        std::vector<T> xs(N);
        typename std::vector<T>::iterator x;
        T val;
        for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
            *x = val;
        return xs;
    }
}

#endif //PROJECT_MATHE_MASTER_MISCELLANEOUS_H
