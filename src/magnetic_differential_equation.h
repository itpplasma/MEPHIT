//
// Created by michi on 03.05.22.
//

#ifndef PHYSIC_MASTER_THESIS_MAGNETIC_DIFFERENTIAL_EQUATION_H
#define PHYSIC_MASTER_THESIS_MAGNETIC_DIFFERENTIAL_EQUATION_H

#include "mfem.hpp"

#include <functional>

namespace mde {
    //assume all functions depend on R and Z, i.e. f = f(R,Z)
    mfem::SparseMatrix
    construct_mde_matrix(mfem::FiniteElementSpace &U, std::function<mfem::Vector(const double &, const double &)> h_t,
                         std::function<double(const double &, const double &)> h_phi, const int n, const int integration_order = -1);

    mfem::Vector construct_mde_rhs(mfem::FiniteElementSpace &U, std::function<mfem::Vector(const double &,
                                                                                           const double &)> h_t,
                                   std::function<double(const double &, const double &)> h_phi,
                                   std::function<double(const double &, const double &)> f_r,
                                   std::function<double(const double &, const double &)> f_i,
                                   const int n, const int integration_order = -1);

    namespace detail{
        mfem::SparseMatrix construct_mde_stiffness(mfem::FiniteElementSpace &U,
        const std::function<mfem::Vector(const double &,
                                         const double &)> &h_t
        , mfem::IntegrationRule* ir);
        mfem::SparseMatrix construct_mde_mass_matrix(mfem::FiniteElementSpace &U,
                                                                  const std::function<double(const double &, const double &)> &h_phi,
                                                                  const int n, const mfem::IntegrationRule *ir);
        mfem::SparseMatrix construct_mde_mixed_term(mfem::FiniteElementSpace &U,
                                                                 const std::function<mfem::Vector(const double &, const double &)> &h_t,
                                                                 const std::function<double(const double &, const double &)> &h_phi,
                                                                 const int n, const mfem::IntegrationRule *ir);
    }

};


#endif //PHYSIC_MASTER_THESIS_MAGNETIC_DIFFERENTIAL_EQUATION_H
