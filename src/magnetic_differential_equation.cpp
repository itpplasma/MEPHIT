//
// Created by michi on 03.05.22.
//

#include "magnetic_differential_equation.h"

#include "miscellaneous.h"


mfem::SparseMatrix mde::detail::construct_mde_stiffness(mfem::FiniteElementSpace &U,
                                                const std::function<mfem::Vector(const double &,
                                                                                 const double &)> &h_t
                                                                                 , mfem::IntegrationRule* ir) {
    //construct stiffness matrix with spacial anisotropy matrix
    std::function<void(const mfem::Vector &, mfem::DenseMatrix &)> H_f = [&h_t](const mfem::Vector & p, mfem::DenseMatrix & Mat){
        assert(p.Size() == 2);
        const double R = p(0);
        const double Z = p(1);
        mfem::Vector h = h_t(R,Z);
        assert(h.Size() == 2);
        //anisotropy matrix defined by M = h_t \cdot \h_t^T
        Mat = mfem::DenseMatrix(2,2);
        Mat(0,0) = h(0)*h(0);
        Mat(0,1) = h(0)*h(1);
        Mat(1,0) = h(1)*h(0);
        Mat(1,1) = h(1)*h(1);
    };
    mfem::MatrixFunctionCoefficient H(2, H_f);
    mfem::BilinearForm k_form(&U);
    k_form.AddDomainIntegrator(new mfem::DiffusionIntegrator(H, ir));
    k_form.Assemble();
    mfem::SparseMatrix K(k_form.SpMat());
    return K;
}

mfem::SparseMatrix mde::detail::construct_mde_mass_matrix(mfem::FiniteElementSpace &U,
                                                  const std::function<double(const double &, const double &)> &h_phi,
                                                  const int n, const mfem::IntegrationRule *ir) {
    std::function<double(const mfem::Vector &)> h_phi_squared_fun= [&h_phi, n](const mfem::Vector& p){
        assert(p.Size() == 2);
        const double R = p(0);
        const double Z = p(1);
        const double f_val = h_phi(R,Z);
        return f_val*f_val*n*n;
    };
    mfem::FunctionCoefficient h_phi_squared_fc(h_phi_squared_fun);
    mfem::BilinearForm mass_form(&U);
    mass_form.AddDomainIntegrator(new mfem::MassIntegrator(h_phi_squared_fc, ir));
    mass_form.Assemble();
    mfem::SparseMatrix M(mass_form.SpMat());
    return M;
}


mfem::SparseMatrix mde::detail::construct_mde_mixed_term(mfem::FiniteElementSpace &U,
                                                 const std::function<mfem::Vector(const double &, const double &)> &h_t,
                                                 const std::function<double(const double &, const double &)> &h_phi,
                                                 const int n, const mfem::IntegrationRule *ir) {
    std::function<void(const mfem::Vector &, mfem::Vector &)> h_phi_h_t_fun = [&h_phi, h_t, n](mfem::Vector const& p, mfem::Vector& res){
        assert(p.Size() == 2);
        const double R = p(0);
        const double Z = p(1);
        mfem::Vector h = h_t(R,Z);
        assert(h.Size() == 2);
        const double mult = h_phi(R,Z)*n;
        h *= -mult; //- due to definition of mfem bilinear form
        res = h;
    };
    mfem::VectorFunctionCoefficient h_phi_h_t(2,h_phi_h_t_fun);
    mfem::MixedBilinearForm g_form(&U, &U);
    mfem::MixedScalarWeakDivergenceIntegrator * integrator = new mfem::MixedScalarWeakDivergenceIntegrator(h_phi_h_t);
    if(ir){
        integrator->SetIntRule(ir);
    }
    g_form.AddDomainIntegrator(integrator);
    g_form.Assemble();
    mfem::SparseMatrix G(g_form.SpMat());
    return G;
}


mfem::SparseMatrix
mde::construct_mde_matrix(mfem::FiniteElementSpace &U, std::function<mfem::Vector(const double &, const double &)> h_t,
                          std::function<double(const double &, const double &)> h_phi, const int n, const int integration_order) {

    mfem::IntegrationRule* ir = nullptr; //use default integration rule
    if(integration_order > 0){
        ir = new mfem::IntegrationRule(mfem::IntRules.Get(mfem::Geometry::Type::TRIANGLE, integration_order));
    }
    mfem::SparseMatrix K = detail::construct_mde_stiffness(U, h_t, ir);
    mfem::SparseMatrix M = detail::construct_mde_mass_matrix(U, h_phi, n, ir);

    K += M;
    K.Finalize();

    mfem::SparseMatrix G = detail::construct_mde_mixed_term(U, h_t, h_phi, n, ir);
    mfem::SparseMatrix G_final = G;
    G_final.Finalize(); //mfem only supports finalized matrices to be transposed
    mfem::SparseMatrix* Gt = mfem::Transpose(G_final);
    G *= -1.;
    G += *Gt;
    delete Gt;

    mfem::SparseMatrix G_minus = G;
    G_minus *= -1;
    G.Finalize();
    G_minus.Finalize();

    //construct block operator
    mfem::Array<int> block_offsets(3);
    block_offsets[0] = 0;
    block_offsets[1] = U.GetVSize();
    block_offsets[2] = U.GetVSize();
    block_offsets.PartialSum();

    mfem::BlockMatrix Sys(block_offsets);
    Sys.SetBlock(0,0, &K);
    Sys.SetBlock(0,1, &G);
    Sys.SetBlock(1,0, &G_minus);
    Sys.SetBlock(1,1, &K);
    Sys.Finalize();

    mfem::SparseMatrix S = utility::mfem_helper::block_to_sparse_matrix(Sys);
    S.Finalize();

    delete ir;

    return S;
}

mfem::Vector
mde::construct_mde_rhs(mfem::FiniteElementSpace &U, std::function<mfem::Vector(const double &, const double &)> h_t,
                       std::function<double(const double &, const double &)> h_phi,
                       std::function<double(const double &, const double &)> f_r,
                       std::function<double(const double &, const double &)> f_i,
                       const int n,
                       const int integration_order) {
    mfem::IntegrationRule* ir = nullptr; //use default integration rule
    if(integration_order > 0){
        ir = new mfem::IntegrationRule(mfem::IntRules.Get(mfem::Geometry::Type::TRIANGLE, integration_order));
    }

    //construct functions for linear form integrators
    std::function<void(const mfem::Vector &, mfem::Vector &)> fr_1_fun = [&h_t, &f_r](mfem::Vector const& p, mfem::Vector& res){
        assert(p.Size() == 2);
        const double R = p(0);
        const double Z = p(1);
        mfem::Vector h = h_t(R,Z);
        h *= f_r(R,Z);
        res = h;
    };
    std::function<void(const mfem::Vector &, mfem::Vector &)> fi_1_fun = [&h_t, &f_i](mfem::Vector const& p, mfem::Vector& res){
        assert(p.Size() == 2);
        const double R = p(0);
        const double Z = p(1);
        mfem::Vector h = h_t(R,Z);
        h *= f_i(R,Z);
        res = h;
    };
    std::function<double(const mfem::Vector &)> fr_2_fun= [&h_phi, &f_r, n](const mfem::Vector& p){
        assert(p.Size() == 2);
        const double R = p(0);
        const double Z = p(1);
        return -f_r(R,Z)*h_phi(R,Z)*n;
    };
    std::function<double(const mfem::Vector &)> fi_2_fun= [&h_phi, &f_i, n](const mfem::Vector& p){
        assert(p.Size() == 2);
        const double R = p(0);
        const double Z = p(1);
        return f_i(R,Z)*h_phi(R,Z)*n;
    };

    mfem::VectorFunctionCoefficient fr_1(2,fr_1_fun);
    mfem::FunctionCoefficient fr_2(fr_2_fun);
    mfem::VectorFunctionCoefficient fi_1(2,fi_1_fun);
    mfem::FunctionCoefficient fi_2(fi_2_fun);


    //assemble linear forms
    mfem::LinearForm first_row(&U);
    first_row.AddDomainIntegrator(new mfem::DomainLFIntegrator(fi_2, ir));
    mfem::DomainLFGradIntegrator* integrator1 = new mfem::DomainLFGradIntegrator(fr_1);
    if(ir) {
        integrator1->SetIntRule(ir);
    }
    first_row.AddDomainIntegrator(integrator1);
    first_row.Assemble();
    mfem::LinearForm second_row(&U);
    second_row.AddDomainIntegrator(new mfem::DomainLFIntegrator(fr_2, ir));
    mfem::DomainLFGradIntegrator* integrator2 = new mfem::DomainLFGradIntegrator(fi_1);
    if(ir){
        integrator2->SetIntRule(ir);
    }
    second_row.AddDomainIntegrator(integrator2);
    second_row.Assemble();

    //put all to rhs
    mfem::Vector rhs(2*U.GetVSize());
    for(int i=0; i<U.GetVSize(); ++i){
        rhs(i) = first_row.Elem(i);
    }
    for(int i=0; i<U.GetVSize(); ++i){
        rhs(i+U.GetVSize()) = second_row.Elem(i);
    }

    delete ir;

    return rhs;
}




