#include "optimization.hpp"
#include "linalg.hpp"
#include <cfloat>


Matrix<float> linprog_can(Matrix<float> f, Matrix<float> A, Matrix<float> b)
{
    std::size_t n_variables = f.rowNb();
    std::size_t n_constraints = b.rowNb();

    Matrix<float> new_f = vertcat(f, zeros<float>(n_constraints,1));
    Matrix<float> new_A = horzcat(A, id<float>(n_constraints));

    return linprog_std(new_f, new_A, b).getRows(0,n_variables);
}

Matrix<float> linprog_std(Matrix<float> f, Matrix<float> A, Matrix<float> b)
{
    std::size_t n_variables = f.rowNb();
    std::size_t n_constraints = A.rowNb();

    Matrix<float> tab(n_constraints, n_constraints+n_variables);
    Matrix<float> Z(1, n_constraints+n_variables);
    Matrix<float> C(n_constraints+1, 1);
    Matrix<float> K = C;

    Matrix<float> new_A = A;
    Matrix<float> new_b = b;
    for(std::size_t i=0;i<n_constraints;++i)
    {
        if(b(i,0)<0.0f)
        {
            new_b(i,0) = -b(i, 0);
            new_A.setRow(i, -A.getRow(i));
        }
    }

    /*Z.setCols(0, -transpose(f));
    tab.setCols(0, new_A);
    tab.setCols(n_variables, id<float>(n_constraints));
    C.setRows(0, new_b);*/
    Matrix<float> Cb = ones<float>(n_constraints, 1)*1e4;
    Z.setCols(0, -transpose(f)+dot(transpose(Cb), new_A));
    tab.setCols(0, new_A);
    tab.setCols(n_variables, id<float>(n_constraints));
    C.setRows(0, new_b);
    C(n_constraints, 0) = -dot(transpose(Cb), new_b)(0,0);

    Matrix<float> X_final = zeros<float>(n_variables, 1);
    Matrix<float> S_final = new_b;
    Matrix<float> outgoing_variable_idx = arange<float>(n_variables, n_variables+n_constraints);

    while(count(Z>0, true)>0)
    {
        std::size_t Cp_idx = argmax(Z)(0,1);
        Matrix<float> Cp = tab.getCol(Cp_idx);
        K = C.getRows(0, n_constraints)*(1.0f/Cp);
        replace_if(K, K<=0, FLT_MAX);
        std::size_t Rp_idx = argmin(K)(0,0);
        float pivot = tab(Rp_idx, Cp_idx);
        tab.setRow(Rp_idx, tab.getRow(Rp_idx)/pivot);
        C(Rp_idx,0)/=pivot;

        for(std::size_t i=0;i<Rp_idx;++i)
        {
            float alpha = tab(i, Cp_idx);
            tab.setRow(i, tab.getRow(i)-alpha*tab.getRow(Rp_idx));
            C(i,0)-=alpha*C(Rp_idx,0);
        }
        for(std::size_t i=Rp_idx+1;i<n_constraints;++i)
        {
            float alpha = tab(i, Cp_idx);
            tab.setRow(i, tab.getRow(i)-alpha*tab.getRow(Rp_idx));
            C(i,0)-=alpha*C(Rp_idx,0);
        }
        float alpha = Z(0, Cp_idx);
        Z-=alpha*tab.getRow(Rp_idx);
        C(n_constraints,0)-=alpha*C(Rp_idx, 0);

        outgoing_variable_idx(Rp_idx, 0) = Cp_idx;

    }

    for(std::size_t i=0;i<n_constraints;++i)
    {
        if(outgoing_variable_idx(i, 0)<n_variables)
            X_final(outgoing_variable_idx(i, 0), 0) = C(i, 0);
        else
            S_final(outgoing_variable_idx(i, 0)-n_variables, 0) = C(i, 0);
    }
    /*std::cout<<outgoing_variable_idx<<std::endl;
    std::cout<<C<<std::endl;
    std::cout<<Z<<std::endl;*/
    return X_final;
}

Matrix<float> linprog(Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq)
{
    Matrix<float> new_A = vertcat(A, vertcat(Aeq, -Aeq));
    Matrix<float> new_b = vertcat(b, vertcat(beq, -beq));

    return linprog_can(f, new_A, new_b);
}


Matrix<float> quadprog(Matrix<float> H, Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq, Matrix<float> X0)
{
    std::size_t N = X0.rowNb();
    std::size_t nb_eq_constraints = Aeq.rowNb();
    std::size_t nb_ineq_constraints = A.rowNb();

    if(X0.colNb()!=1)
    {
        std::cout<<"X0 must be a column vector"<<std::endl;
        return X0;
    }
    else if((H.rowNb() != N) || (H.colNb() != N))
    {
        std::cout<<"H must be a "<<N<<"x"<<N<<" matrix"<<std::endl;
        return X0;
    }
    else if((f.rowNb() != N) || (f.colNb() != 1))
    {
        std::cout<<"f must be a "<<N<<"x"<<"1"<<" matrix"<<std::endl;
        return X0;
    }
    else if(A.colNb()!=N)
    {
        std::cout<<"A must have "<<N<<" columns"<<std::endl;
        return X0;
    }
    else if((A.rowNb()!=b.rowNb()) || (b.colNb()!=1))
    {
        return X0;
    }
    else if(Aeq.colNb()!=N)
    {
        std::cout<<"A must have "<<N<<" columns"<<std::endl;
        return X0;
    }
    else if((Aeq.rowNb()!=beq.rowNb()) || (beq.colNb()!=1))
    {
        return X0;
    }
    else
    {
        Matrix<float> x = X0;
        Matrix<bool> mask = ones<bool>(nb_ineq_constraints, 1);
        Matrix<float> A_t = transpose(A);
        Matrix<float> Aeq_t = transpose(Aeq);
        std::size_t nb_active_ineq_constraints = nb_ineq_constraints;

        while(true)
        {
            Matrix<float> TKL_A(N+nb_eq_constraints+nb_active_ineq_constraints);
            Matrix<float> TKL_B(N+nb_eq_constraints+nb_active_ineq_constraints, 1);

            TKL_A.setSubmat(0, 0, Aeq);
            if(nb_active_ineq_constraints>0)
                TKL_A.setSubmat(nb_eq_constraints, 0, compress(mask, A, 1));
            TKL_A.setSubmat(nb_eq_constraints+nb_active_ineq_constraints, 0, H);
            TKL_A.setSubmat(nb_eq_constraints+nb_active_ineq_constraints, N, Aeq_t);
            if(nb_active_ineq_constraints>0)
                TKL_A.setSubmat(nb_eq_constraints+nb_active_ineq_constraints, N+nb_eq_constraints, compress(mask, A_t, 2));

            TKL_B.setRows(0,beq-dot(Aeq, x));
            if(nb_active_ineq_constraints>0)
                TKL_B.setRows(nb_eq_constraints,compress(mask, b, 1)-dot(compress(mask, A, 1), x));
            TKL_B.setRows(nb_eq_constraints+nb_active_ineq_constraints,-f-dot(H, x));

            Matrix<float> TKL_X = dot(pinv(TKL_A), TKL_B);
            Matrix<float> d = TKL_X.getRows(0, N);
            Matrix<float> lambda = TKL_X.getRows(N, N+nb_eq_constraints);
            Matrix<float> mu;
            if(nb_active_ineq_constraints>0)
                mu = TKL_X.getRows(N+nb_eq_constraints, N+nb_eq_constraints+nb_active_ineq_constraints);

            std::cout<<TKL_X<<std::endl;
            std::cout<<x<<std::endl;
            if(norm(d)<10e-9)
            {
                    if(nb_active_ineq_constraints==0 || count(mu>=0,true)==nb_active_ineq_constraints)
                        return x;
                    else
                    {
                        Matrix<std::size_t> argminj = argmin(mu);//FAUX
                        --nb_active_ineq_constraints;
                        mask(argminj(0, 0), 1) = false;
                    }
            }
            else
            {
                if(count(dot(A, x+d)>=b, true)==nb_ineq_constraints)
                    x+=d;
                else
                {
                    float t = 1.0f;
                    std::cout<<"not implemented yet"<<std::endl;
                    x+=t*d;
                }
            }
        }
    }
}
