#include "optimization.hpp"
#include "linalg.hpp"
#include <cfloat>
#include <deque>


Matrix<float> steepest_descent(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0)
{
    Matrix<float> X = X0;
    Matrix<float> grad = ones<float>(X.rowNb(),X.colNb());
    std::size_t cpt = 0;
    while(norm(grad)>10e-10 && cpt<100)
    {
        auto evalfunc = (*ptr)(X);
        //std::cout<<"X "<<X<<std::endl;
        float func = evalfunc.first;
        std::cout<<"err "<<func<<std::endl;
        grad = evalfunc.second;
        //std::cout<<"grad "<<grad<<std::endl;
        float alpha = WolfeLineSearch(ptr, X, func, grad, -grad);
        //std::cout<<"alpha "<<alpha<<std::endl;
        X-=alpha*grad;
        ++cpt;
    }
    return X;
}

Matrix<float> bfgs(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0)
{

    Matrix<float> X = X0;
    Matrix<float> H = id<float>(X.size());
    std::size_t cpt = 0;
    auto evalfunc = (*ptr)(X);
    float func_old = evalfunc.first;
    Matrix<float> grad_old = evalfunc.second;
    Matrix<float> grad = grad_old;

    while(norm(grad)>10e-10 && cpt<100)
    {
        Matrix<float> p = -dot(H,grad_old);
        float alpha = WolfeLineSearch(ptr, X, func_old, grad_old, p);
        X+=alpha*p;

        auto evalfunc = (*ptr)(X);
        float func = evalfunc.first;
        grad = evalfunc.second;
        std::cout<<"err "<<func<<std::endl;
        Matrix<float> s = alpha*p;
        Matrix<float> y = grad-grad_old;
        float sy = sum(y*s);
        if(std::abs(sy)<10e-9)
            break;
        Matrix<float> temp = id<float>(X.size())-dot(s,transpose(y))/sy;

        if(cpt==0)
            H = sum(y*s)/sum(y*y)*id<float>(X.size());

        H = dot(temp,dot(H,transpose(temp)))+dot(s,transpose(s))/sy;

        func_old = func;
        grad_old = grad;

        ++cpt;
    }
    return X;
}

Matrix<float> lbfgs(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0)
{
    std::cout<<"hello"<<std::endl;
    int m = 5;
    std::deque<Matrix<float> > s_queue, y_queue;

    Matrix<float> X = X0;
    std::size_t cpt = 0;

    auto evalfunc = (*ptr)(X);
    float func_old = evalfunc.first;
    Matrix<float> grad_old = evalfunc.second;
    Matrix<float> grad = grad_old;

    while(norm(grad)>10e-10 && cpt<100)
    {
        std::cout<<cpt<<std::endl;
        Matrix<float> p = grad_old;
        if(cpt>0)
        {
            Matrix<float> cur_y, cur_s;
            std::deque<float> alpha_queue,rho_queue;
            for(std::size_t i=0;i<s_queue.size();++i)
            {
                cur_y = y_queue[i];
                cur_s = s_queue[i];
                float rho = sum(cur_y*cur_s);
                std::cout<<"rho"<<rho<<std::endl;
                rho_queue.push_back(rho);
                float alpha = sum(cur_s*p)/rho;
                alpha_queue.push_back(alpha);
                p-=alpha*cur_y;
            }
            p *= rho_queue.front()/sum(y_queue[0]*y_queue[0]);
            std::cout<<p<<std::endl;
            for(std::size_t i=s_queue.size();i>0;--i)
            {
                float beta = sum(p*y_queue[i-1])/rho_queue[i-1];
                p+=s_queue[i-1]*(alpha_queue[i-1]-beta);
            }
            alpha_queue.clear();
            rho_queue.clear();
        }
        p = -p;
        std::cout<<p<<std::endl;
        float alpha = WolfeLineSearch(ptr, X, func_old, grad_old, p);
        X+=alpha*p;

        auto evalfunc = (*ptr)(X);
        float func = evalfunc.first;
        grad = evalfunc.second;
        s_queue.push_front(alpha*p);
        y_queue.push_front(grad-grad_old);

        func_old = func;
        grad_old = grad;
        if(s_queue.size()>m)
        {
            s_queue.pop_back();
            y_queue.pop_back();
        }
        ++cpt;
    }
    return X;
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

Matrix<float> linprog_can(Matrix<float> f, Matrix<float> A, Matrix<float> b)
{
    std::size_t n_variables = f.rowNb();
    std::size_t n_constraints = b.rowNb();

    Matrix<float> new_f = vertcat(f, zeros<float>(n_constraints,1));
    Matrix<float> new_A = horzcat(A, id<float>(n_constraints));

    return linprog_std(new_f, new_A, b).getRows(0,n_variables);
}

Matrix<float> linprog(Matrix<float> f, Matrix<float> A, Matrix<float> b, Matrix<float> Aeq, Matrix<float> beq)
{
    std::size_t n_variables = f.rowNb();
    std::size_t n_constraints = b.rowNb();

    Matrix<float> new_f = vertcat(f, zeros<float>(n_constraints,1));
    Matrix<float> new_Aeq = vertcat(horzcat(Aeq,zeros<float>(n_constraints)), horzcat(A,id<float>(n_constraints)));
    Matrix<float> new_beq = vertcat(beq, b);

    return linprog_std(new_f, new_Aeq, new_beq).getRows(0,n_variables);
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

// predefined function
std::pair<float, Matrix<float> > rosenbrock(const Matrix<float>& X)
{
    float x = X(0,0);
    float y = X(1,0);
    float f = (1.0f-x)*(1.0f-x)+100.0f*(y-x*x)*(y-x*x);
    Matrix<float> g(2,1);
    g(0,0) = -2.0f*(1.0f-x)-400.0f*x*(y-x*x);
    g(1,0) = 200.0f*(y-x*x);
    return std::make_pair(f, g);
}

std::pair<float, Matrix<float> > himmelblau(const Matrix<float>& X)
{
    float x = X(0,0);
    float y = X(1,0);
    float f = std::pow(x*x+y-11.0f, 2.0f)+std::pow(x+y*y-7.0f, 2.0f);
    Matrix<float> g(2,1);
    g(0,0) = 2.0f*(2.0f*x*(x*x+y-11.0f)+(x+y*y-7.0f));
    g(1,0) = 2.0f*((x*x+y-11.0f)+2.0f*y*(x+y*y-7.0f));
    return std::make_pair(f, g);
}

std::pair<float, Matrix<float> > rastrigin(const Matrix<float>& X)
{
    float x = X(0,0);
    float y = X(1,0);
    float f = 10.0f + x*x-10.0f*cos(2.0f*PI*x) + y*y-10.0f*cos(2.0f*PI*y);
    Matrix<float> g(2,1);
    g(0,0) = 2.0f*x + 10.0f*2.0f*PI*sin(2.0f*PI*x);
    g(1,0) = 2.0f*y + 10.0f*2.0f*PI*sin(2.0f*PI*y);
    return std::make_pair(f, g);
}

float ArmijoLineSearch(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p)
{
    float alpha = 1.0f;
    float coeff = 0.2f;
    float c1 = 0.1f;
    float dir = sum(p*g0);
    while(true)
    {
        float f_x_ap = ((*ptr)(X0+alpha*p)).first;
        if(f_x_ap>err+c1*alpha*dir)
            alpha*=coeff;
        else
            return alpha;
    }
}

float WolfeLineSearch(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X, float err, const Matrix<float>& g0, const Matrix<float>& p)
{
    float alpha_min = 0.0f;
    float alpha_max = FLT_MAX;
    float t = 1.0f;
    float c1 = 0.1f;
    float c2 = 0.9f;
    float dir = sum(p*g0);
    while(true)
    {
        auto f_x_ap = (*ptr)(X+t*p);
        if(f_x_ap.first > err+c1*t*dir)
        {
            alpha_max = t;
            t = 0.5f*(alpha_min+alpha_max);
        }
        else if(sum(p*f_x_ap.second) < c2*dir)
        {
            alpha_min = t;
            if(alpha_max==FLT_MAX)
                t = 2.0f*alpha_min;
            else
                t = 0.5f*(alpha_min+alpha_max);
        }
        else
        {
            return t;
        }
    }
}

Matrix<float> conjugate_gradient(const Matrix<float>& A, const Matrix<float>& b, const Matrix<float>& X)
{
    Matrix<float> X0 = X;
    Matrix<float> r = dot(A, X0)-b;
    Matrix<float> p = -r;
    std::size_t k = 0;
    std::size_t ndims = X0.size();
    float r2 = sum(r*r);
    while(k<ndims || r2>1e-10)
    {
        Matrix<float> Ap = dot(A, p);
        float alpha = r2/sum(p*Ap);
        X0 += alpha*p;
        r += alpha*Ap;
        float beta = sum(r*r)/r2;
        r2 = sum(r*r);
        p = -r+beta*p;
        ++k;
    }
    return X0;
}

float zoom(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p, float alphal, float alphah)
{
    float c1 = 1e-4;
    float c2 = 0.5f;

    auto f_x0 = (*ptr)(X0);
    float fx0 = f_x0.first;
    Matrix<float> gx0 = f_x0.second;
    float slope0 = sum(p*gx0);

    auto f_x_lo = (*ptr)(X0+alphal*p);
    float fxlo = f_x_lo.first;

    std::size_t i=0;
    while(true)
    {
        float alpha = (alphal+alphah)/2.0f;

        Matrix<float> Xc = X0+alpha*p;
        auto f_x_ap = (*ptr)(Xc);
        float fxc = f_x_ap.first;
        Matrix<float> gxc = f_x_ap.second;

        if((fxc>fx0+c1*alpha*slope0)||(fxc>=fxlo))
        {
            alphah = alpha;
        }
        else
        {
            float slopexc = sum(p*gxc);
            if(std::abs(slopexc)<=-c2*slope0)
            {
                return alpha;
            }
            if(slopexc*(alphah-alphal)>=0)
            {
                alphah = alphal;
            }
            alphal = alpha;
            fxlo = fxc;
        }

        ++i;
    }
}



float StrongWolfeLineSearch(std::pair<float, Matrix<float> > (*ptr)(const Matrix<float>&), const Matrix<float>& X0, float err, const Matrix<float>& g0, const Matrix<float>& p)
{
    float c1 = 1e-4;
    float c2 = 0.5f;

    float alphamin = 0;
    float alphamax = FLT_MAX;
    float alpha_old = 0;
    float alpha = 1;

    auto f_x0 = (*ptr)(X0);
    float fx0 = f_x0.first;
    Matrix<float> gx0 = f_x0.second;
    float slope0 = sum(p*gx0);

    std::size_t i=0;

    while(true)
    {

        Matrix<float> Xc = X0+alpha*p;
        auto f_x_ap = (*ptr)(Xc);
        float fxc = f_x_ap.first;
        Matrix<float> gxc = f_x_ap.second;
        float slopexc = sum(p*gxc);

        if((fxc>fx0+c1*alpha*slope0)||(fxc>=fx0 && i>0))
        {
            return zoom(ptr, X0,err, g0, p, alpha_old, alpha);
        }
        if(std::abs(slopexc)<=-c2*slope0)
        {
            return alpha;
        }
        if(slopexc >= 0.0f)
        {
            return zoom(ptr, X0,err, g0, p, alpha_old, alpha);
            return alpha;
        }

        alpha_old = alpha;
        alpha = 2*alpha;

        ++i;
    }
}
