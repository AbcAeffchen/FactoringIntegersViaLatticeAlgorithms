
#include "NewEnum.h"

NewEnumStage::NewEnumStage(const RR& y_t, const RR& c_t, const RR& c_tp1, const Vec<RR>& u, unsigned long t, Vec<RR> target, RR v_adjust)
{
    this->y_t = y_t;
    this->c_t = c_t;
    this->c_tp1 = c_tp1;
    this->u = conv<Vec<double> >(u);
    this->t = t;
    this->target = target;
    this->v_adjust = v_adjust;
}

Vec<RR> NewEnumStage::get_u()
{
    return conv<Vec<RR> >(this->u);
}

void NewEnum::precomputeLevelProbabilities()
{
    this->level_probabilitys.SetLength(this->s_max);

    for(long level = 10; level <= this->s_max; level++)
    {
        this->level_probabilitys(level) = power2_RR(level * (-1));
    }
}

RR NewEnum::closest_RR (RR x)
{
    return ceil(x-0.5);
}

RR NewEnum::next(RR u, RR y)
{
    long side1 = sign(this->closest_RR(y) - y);
    if(side1 == 0)
        side1 = 1;

    long side2 = sign(u - y);
    if(side2 == 0)
        side2 = 1;

    RR next = 2 * this->closest_RR(y) - u;

    if(side1 == side2 || this->closest_RR(y) == u)
        return next - side2;
    else
        return next;
}

void NewEnum::clearL()
{
    vector<queue<NewEnumStage> > empty;
    swap(this->L,empty);
}

void NewEnum::perform()
{
    // Reset list of delayed stages
    this->current_s = 10;
    this->clearL();
    this->L = vector<queue<NewEnumStage> > (this->s_max - 9);          // setting the number of queues of stages

    Vec<RR> u;                              // coordinates of a close vector
    u.SetLength(this->n);

    u(this->n) = this->closest_RR(this->tau(this->n));

    NewEnumStage start = NewEnumStage(this->tau(this->n), conv<RR>(0), conv<RR>(0),u, this->n, this->tau, conv<RR>(1));
    this->perform(start, true);
}

void NewEnum::perform(NewEnumStage& current_stage, bool perform_delayed_stages)
{
    this->decrease_max_distance = true;

    ZZ u_rel, v_rel;

    long max_t = current_stage.t;
    long t = current_stage.t;                   // projection to the n-t-1 last coordinates
    Vec<RR> c,                          // Length in the projection
            u,                          // coordinates of a close vector
            y;                          // uncovered parts of the target vector
    c.SetLength(t+1);
    c(t) = current_stage.c_t;
    c(t+1) = current_stage.c_tp1;
    u = current_stage.get_u();
    y.SetLength(t);
    y(t) = current_stage.y_t;

    RR reduce_ratio;

    RR rho;
    RR beta;
    long level;

    bool success;

    // perform stages with s = current_s
    while(t <= max_t)
    {
        c(t) = c(t+1) + power(abs(u(t) - y(t)), 2) * this->R_ii_squared(t);

        if(c(t) >= this->A_curr)
        {
            // 2.1
            goto cleanup;
        }

        if(t == 1)
        {
            success = this->getEquation(u,c(1));

            if (c(1) < this->A_curr)
            {
                reduce_ratio = c(1) / this->A_curr; // 0 < reduce_factor < 1
                if(this->decrease_max_distance || reduce_ratio < this->min_reduce_ratio)
                {
                    this->A_curr = c(1);        // reduce max distance
//                        if(this->current_s == 10 && reduce_ratio < this->min_restart_ratio)      // min_restart_factor = 0 -> never restart, min_restart_factor = 1 -> always restart
                    if(this->current_s == 10)      // min_restart_factor = 0 -> never restart, min_restart_factor = 1 -> always restart
                    {
                        this->perform();        // neustart
                        return;
                    }
                }
                else
                {
                    // do nothing... works fine
                }
            }

            if(success) // check for equation
            {
                // do statistics here

                this->decrease_max_distance = false;

            }
            // 2.1
            goto cleanup;
        }


        // compute the probability beta
        rho = SqrRoot(power(this->A_curr - c(t), t-1));
        beta = this->V(t-1) * rho;
        if(beta >= this->level_probabilitys(this->current_s) * conv<RR>(t))
        {
            t--;

            clear(y(t));        // y(t) = 0
            for(long i = t + 1; i <= this->n; i++)    // 1/r_tt * sum_{i=t+1}^n (\tau_i - u_i) r_ti
                y(t) += (current_stage.target(i) - u(i)) * this->mu(i,t);
            y(t) += current_stage.target(t);

            u(t) = this->closest_RR(y(t));
            continue;
        }
        else
        {
            level = conv<long>(ceil((- log(beta) + log(t))/log(2)));
            if(level <= this->s_max)
            {
                this->L[level - 10].emplace(y(t),c(t),c(t+1),u,t,current_stage.target,current_stage.v_adjust);      // save the stage for later
            }
            // 2.1
            goto cleanup;
        }

        // 2.1
        cleanup:
        t++;
        if(t > max_t)
            break;
        u(t) = this->next(u(t), y(t));
    }

    // 3
    // perform the stages with 10 < s <= s_max
    long delayed_before = 0;
    long delayed_now = 0;
    if(perform_delayed_stages)
    {
        this->delayedStagesCounter[0] = 0;
        while(this->current_s <= this->s_max)
        {
            // Output of delayed stages status

            cout << "#L[" << (this->current_s - 10) << "] = " << this->L[this->current_s - 10].size();
            this->delayedStagesCounter[this->current_s - 10] = this->L[this->current_s - 10].size();
            this->delayedStagesCounter[0] += this->L[this->current_s - 11].size();
            delayed_now = 0;
            for(int i = 0; i < this->s_max - 10; i++)
                delayed_now += this->L[i].size();

            cout << "  |  " << delayed_now << "  |  " << (delayed_now - delayed_before) << endl;
            delayed_before = delayed_now;

            // perform delayed stages
            while(!this->L[this->current_s - 10].empty())
            {
                this->perform(this->L[this->current_s - 10].front(), false);
                this->L[this->current_s - 10].pop();
            }
            this->current_s++;
        }
    }
}

void NewEnum::precomputeVolumes(long dim)
{
    RR pi = ComputePi_RR();
    this->V.SetLength(dim);
    Vec<RR> R_diag_prod;
    R_diag_prod.SetLength(dim);

    if(dim >= 1)
        this->V(1) = conv<RR>(2);
    if(dim >= 2)
        this->V(2) = pi;

    for(long i = 3; i <= dim; i++)
    {
        this->V(i) = this->V(i-2) * pi / (i/2.0);
    }

    R_diag_prod(1) = SqrRoot(this->R_ii_squared(1));
    this->V(1) /= R_diag_prod(1);
    for(long i = 2; i <= dim; i++)
    {
        R_diag_prod(i) = SqrRoot(this->R_ii_squared(i)) * R_diag_prod(i-1);
        this->V(i) /= R_diag_prod(i);
    }

    return;
}

void NewEnum::computeUV(const Vec<RR>& input, ZZ& u, ZZ& v)
{
    Vec<RR> close_vec;

    Vec<RR> temp;
    mul(temp, conv<Mat<RR> >(this->U_scaled), input);
    temp += this->shift;
    mul(close_vec, conv<Mat<RR> >(this->U), temp);

    Vec<long> equation;
    v = conv<ZZ>(1);
    equation.SetLength(close_vec.length() + 1);     // exponents of the first primes and -1

    u = conv<ZZ>(1);

    for(long i = 1; i <= close_vec.length(); i++)
    {
        if(close_vec(i) > 0)
        {
            u *= power_ZZ(this->primes(i), conv<long>(close_vec(i)));
            equation(i) = conv<long>(close_vec(i));
        }
        else
        {
            equation(i) = 0;
            v *= power_ZZ(this->primes(i), conv<long>(-close_vec(i)));
        }
    }
}

bool NewEnum::getEquation(const Vec<RR>& input, RR& c_1)
{
    mul(this->temp_vec, conv<Mat<RR> >(this->U_scaled), input);
    this->temp_vec += this->shift;
    mul(this->close_vec, conv<Mat<RR> >(this->U), this->temp_vec);

    this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
    this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1

    NTL::set(this->u);       // u = 1

    for(long i = 1; i <= this->n; i++)
    {
        if(this->close_vec(i) > 0)
        {
            this->u *= power_ZZ(this->primes(i), conv<long>(this->close_vec(i)));
            this->raw_equation(i) = conv<long>(this->close_vec(i));
        }
        else
        {
            this->raw_equation(i) = 0;
        }
    }

    this->v = conv<RR>(this->u) / conv<RR>(this->N);
    this->d = this->v - this->closest_RR(this->v);
    this->alpha_nm1 = NTL::abs(this->d);
    this->vN = conv<ZZ>(this->closest_RR(this->v)) * this->N;
    NTL::clear(this->h_n); NTL::set(this->h_nm1); NTL::clear(this->h_nm2);
    NTL::set(this->k_n); NTL::clear(this->k_nm1); NTL::set(this->k_nm2);
    NTL::clear(this->a_nm1);
    long sign = NTL::sign(this->d),
         equation_counter = 0;

    do
    {
        this->nextContinuedFraction(this->h_n, this->k_n, this->h_nm1, this->k_nm1, this->h_nm2, this->k_nm2, this->a_nm1, this->alpha_nm1);

        this->equation = this->raw_equation;
        this->left_side = this->u * this->k_n;
        this->right_side = abs(this->left_side - this->vN * this->k_n - sign * this->h_n * this->N);

        if(this->isSmooth(this->equation, this->k_n, this->left_side, this->right_side))
        {
            equation_counter++;
            this->ride_side_factor = conv<ZZ>(this->closest_RR(conv<RR>(left_side) / conv<RR>(this->N)));
            this->equations.emplace_back(this->equation,this->ride_side_factor,this->current_s,conv<double>(c_1 / this->theoreticalMaxDistance),this->round, this->timer.step());
        }
    }
    while(this->k_n < this->threshold);

    return equation_counter > 0;
}

NewEnum::NewEnum(Timer& timer, FileOutput& file, Statistics& stats, long round, ZZ N, const Vec<long>& primes, const Mat<ZZ>& basis,
    const Mat<ZZ>& U, const Mat<ZZ>& U_scaled, Vec<RR> target_coordinates, Vec<RR> target_shift,
    int s_max, double min_restart_ratio, double min_reduce_ratio, double start_factor_A )
: timer(timer), file(file), stats(stats), B(basis), U(U), U_scaled(U_scaled), primes(primes)
{
    this->round = round;
    this->N = N;
    this->decrease_max_distance = true;
    this->n = this->primes.length();

    // setting up the getEquation workspace
    this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
    this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1
    this->threshold = power_ZZ(this->primes(this->n),3),

    ComputeGS(transpose(this->B),this->mu, this->R_ii_squared);

    this->tau = target_coordinates;
    this->shift = target_shift;

    // Setting the maximum pruning level
    if(s_max > 10)
        this->s_max = s_max;
    else
        this->s_max = 10;

    this->delayedStagesCounter.resize(this->s_max-9,0);     // 0 is used for total counter

    if(min_restart_ratio < 0)
        this->min_restart_ratio = 0;        // never restart
    else
        this->min_restart_ratio = min_restart_ratio;

    if(min_reduce_ratio < 0)
        this->min_reduce_ratio = 0;        // never reduce
    else
        this->min_reduce_ratio = min_reduce_ratio;

    // setting the decreasing behavior
    this->decrease_max_distance = true;

    // precompute (Ballvolume...)
    this->precomputeVolumes(this->R_ii_squared.length());
    this->precomputeLevelProbabilities();

    // start algorithm with a start parameter A
    this->A_curr = 0;
    for(long i = 1; i <= this->R_ii_squared.length(); i++)
    {
        this->A_curr += this->R_ii_squared(i);
    }

    this->A_curr *= 0.25;
    this->theoreticalMaxDistance = this->A_curr;
    this->A_curr *= start_factor_A;
    RR heuristicMaxDistance = this->A_curr;

    this->perform();

    this->stats.updateRoundStats(this->delayedStagesCounter[0] > 0, this->equations.size() > 0);
    this->stats.updateDistanceStats(this->theoreticalMaxDistance, heuristicMaxDistance, this->A_curr);
    this->file.statisticsDistances(this->theoreticalMaxDistance, heuristicMaxDistance, this->A_curr);
    this->file.statisticsDelayedStagesOnLevel(this->delayedStagesCounter);
}

list<Equation> NewEnum::getEquations()
{
    return this->equations;
}


void NewEnum::nextContinuedFraction(ZZ &h_n, ZZ &k_n, ZZ &h_nm1, ZZ &k_nm1, ZZ &h_nm2, ZZ &k_nm2, ZZ &a_nm1, RR &alpha_nm1)
{
    h_n = a_nm1 * h_nm1 + h_nm2;
    k_n = a_nm1 * k_nm1 + k_nm2;
    h_nm2 = h_nm1;
    h_nm1 = h_n;
    k_nm2 = k_nm1;
    k_nm1 = k_n;
    alpha_nm1 = inv(alpha_nm1 - floor(alpha_nm1));       // this is actually alpha_n
    a_nm1 = FloorToZZ(alpha_nm1);                        // this is actually a_n
}

bool NewEnum::isSmooth(Vec<long>& equation, ZZ& k_n, ZZ& left_side, ZZ& right_side)
{
    // check smoothness of k_n
    for(long i = 1; i <= this->n; i++)
    {
        while(k_n % this->primes(i) == 0)
        {
            equation(i) += 1;
            k_n /= this->primes(i);
        }
    }

    if(k_n > 1)        // if not smooth
        return false;

    if(right_side < 0)
        equation(this->n + 1) = 1;
    else
        equation(n + 1) = 0;

    ZZ right_side_abs = abs(right_side);

    // check smoothness of the right side
    for(long i = 1; i <= this->n; i++)
    {
        while(right_side_abs % this->primes(i) == 0)
        {
            if(equation(i) > 0)       // cancel
                left_side /= this->primes(i);
            equation(i) -= 1;
            right_side_abs /= this->primes(i);
        }
    }

    return !(right_side_abs > 1 || left_side == 1);     // if not smooth or the equation is 1 = 1
}