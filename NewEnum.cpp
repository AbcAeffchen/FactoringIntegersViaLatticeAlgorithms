
#include "NewEnum.h"

void NewEnumStage::get_u(Vec<RR> &u) const
{
    conv(u,this->u);
}

Vec<double> NewEnum::precomputeLogT(long n)
{
    Vec<double> log_t;
    log_t.SetLength(n);
    for(long t = 1; t <= n; t++)
        log_t(t) = log(t);

    return log_t;
}

RR NewEnum::closest_RR(const RR &x)
{
    return ceil(x-0.5);
}

void NewEnum::closest_RR(RR &out, const RR &x)
{
    ceil(out, x-0.5);
}

void NewEnum::next(RR &out, const RR &u, const RR &y)
{
    RR closest_y;
    this->closest_RR(closest_y,y);

    int side1 = closest_y >= y ? 1 : -1;
    int side2 = u >= y ? 1 : -1;

    out = 2 * closest_y - u;

    if(side1 == side2 || closest_y == u)
        out -= side2;
}

void NewEnum::clearL()
{
    vector<queue<NewEnumStage>> empty = vector<queue<NewEnumStage>>(this->max_level - 10);
    swap(this->L,empty);
}

void NewEnum::run(unsigned long round, const Mat<ZZ> &newBasis, const Mat<ZZ> &new_U_scaled,
                  const Vec<RR> &new_target_coordinates)
{
    this->prepare(round, newBasis, new_U_scaled, new_target_coordinates);

    cout << "Performing: ";
    // Reset list of delayed stages
    this->current_level = 10;
    this->clearL();

    Vec<RR> u;                              // coordinates of a close vector
    u.SetLength(this->n);
    this->closest_RR(u(this->n),this->tau(this->n));

    this->perform(NewEnumStage(this->tau(this->n), conv<RR>(0), conv<RR>(0),u, this->n));

    this->delayedStagesCounter[0] = 0;

    for(this->current_level = 11; this->current_level <= this->max_level; this->current_level++)
    {
        this->delayedStagesCounter[this->current_level - 10] = this->L[this->current_level - 11].size();
        this->delayedStagesCounter[0] += this->L[this->current_level - 11].size();

        // perform delayed stages
        while(!this->L[this->current_level - 11].empty())
        {
            this->perform(this->L[this->current_level - 11].front());
            this->L[this->current_level - 11].pop();
        }

        cout << ".";
    }

    cout << endl;

    // statistics
    this->stats.updateRoundStats(this->delayedStagesCounter[0] > 0, this->equations.size() > 0);
    this->stats.updateDistanceStats(this->theoreticalMaxDistance, this->heuristicMaxDistance, this->A_curr);
    this->file.statisticsDistances(this->theoreticalMaxDistance, this->heuristicMaxDistance, this->A_curr);
    this->file.statisticsDelayedStagesOnLevel(this->delayedStagesCounter);

}

void NewEnum::perform(const NewEnumStage& current_stage)
{
    this->decrease_max_distance = true;

    long max_t = current_stage.t;
    long t = current_stage.t;           // projection to the n-t-1 last coordinates
    Vec<RR> c,                          // Length in the projection
            u,                          // coordinates of a close vector
            y;                          // uncovered parts of the target vector
    c.SetLength(t+1);
    c(t) = current_stage.c_t;
    c(t+1) = current_stage.c_tp1;
    current_stage.get_u(u);
    y.SetLength(t);
    y(t) = current_stage.y_t;

    RR reduce_ratio;

    long level;

    bool success;

    // perform stages with s = current_level
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
            success = this->checkForEquation(u, c(1));

            if(c(1) < this->A_curr)
            {
                if(this->decrease_max_distance || c(1) / this->A_curr < this->min_reduce_ratio)
                {
                    this->A_curr = c(1);        // reduce max distance
                    // todo add stage level recalculation
                    if(this->current_level == 10)
                    {
                        // do something???
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
        if(t < 4 || (level = calculateLevel(t,c(t))) <= this->current_level)
        {
            t--;

            clear(y(t));        // y(t) = 0
            for(long i = t + 1; i <= this->n; i++)    // 1/r_tt * sum_{i=t+1}^n (\tau_i - u_i) r_ti
                y(t) += (this->tau(i) - u(i)) * this->mu(i,t);
            y(t) += this->tau(t);

            this->closest_RR(u(t),y(t));
            continue;
        }

        // the program comes only this far, if level > current_s holds
        if(level <= this->max_level)
        {
            this->L[level - 11].emplace(y(t),c(t),c(t+1),u,t);      // store the stage for later
        }

        // 2.1
        cleanup:
        t++;
        if(t > max_t)
            break;
        this->next(u(t), u(t), y(t));
    }
}

Vec<double> NewEnum::precomputeVolumes(long n)
{
    Vec<RR> V;
    Vec<double> log_V;
    V.SetLength(n);
    log_V.SetLength(n);

    RR pi = ComputePi_RR();
    V.SetLength(n);
    Vec<RR> R_diag_prod;
    R_diag_prod.SetLength(n);

    if(n >= 1)
        V(1) = conv<RR>(2);
    if(n >= 2)
        V(2) = pi;

    for(long i = 3; i <= n; i++)
    {
        V(i) = V(i-2) * pi / (i/2.0);
        log_V(i) = log(conv<double>(V(i)));
    }

    return log_V;
}

void NewEnum::precomputeLogV()
{
    Vec<double> log_R_diag_prod;
    conv(log_R_diag_prod, this->R_ii_squared);

    for(long i = 1; i <= this->n; i++)
    {
        log_R_diag_prod(i) = log(log_R_diag_prod(i)) / 2.0;
    }

    for(long i = 2; i <= this->n; i++)
    {
        log_R_diag_prod(i) += log_R_diag_prod(i - 1);
    }

    this->log_V_minus_log_R_prod = this->log_V;

    for(long i = 1; i <= this->n; i++)
    {
        this->log_V_minus_log_R_prod(i) -= log_R_diag_prod(i);
    }
}

bool NewEnum::checkForEquation(const Vec<RR> &input, RR &c_1)
{
    mul(this->temp_vec, this->U_scaled_RR, input);
    this->temp_vec += this->shift;
    mul(this->close_vec, this->U_RR, this->temp_vec);

    conv(this->close_vec_long,this->close_vec);

    this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
    this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1

    NTL::set(this->u);       // u = 1

    ZZ temp_ZZ;
    RR temp_RR;

    for(long i = 1; i <= this->n; i++)
    {
        if(this->close_vec(i) > 0)
        {
            power(temp_ZZ, this->primes(i), this->close_vec_long(i));
            this->u *= temp_ZZ;
            this->raw_equation(i) = this->close_vec_long(i);
        }
        else
        {
            this->raw_equation(i) = 0;
        }
    }

//    this->v = conv<RR>(this->u) / conv<RR>(this->N);
    conv(this->v,this->u); conv(temp_RR,this->N);
    this->v /= temp_RR;

//    this->d = this->v - this->closest_RR(this->v);
    this->closest_RR(temp_RR,this->v);
    this->d = this->v - temp_RR;

    NTL::abs(this->alpha_nm1,this->d);

//    this->vN = conv<ZZ>(this->closest_RR(this->v)) * this->N;
    this->closest_RR(temp_RR,this->v);
    conv(this->vN,temp_RR);
    this->vN *= this->N;

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
        abs(this->right_side,this->left_side - this->vN * this->k_n - sign * this->h_n * this->N);

        if(this->isSmooth(this->equation, this->k_n, this->left_side, this->right_side))
        {
            equation_counter++;
            this->ride_side_factor = conv<ZZ>(this->closest_RR(conv<RR>(left_side) / conv<RR>(this->N)));
            this->equations.emplace_back(this->equation, this->ride_side_factor, this->current_level, conv<double>(c_1 / this->theoreticalMaxDistance), this->round, this->timer.step());
        }
    }
    while(this->k_n < this->threshold);

    return equation_counter > 0;
}

NewEnum::NewEnum(const FactoringSettings &settings, Timer &timer, FileOutput &file,
                 Statistics &stats, const Vec<long> &primes, const Mat<ZZ> &U, const Vec<RR> &target_shift)
        : settings(settings), timer(timer), file(file), stats(stats), N(settings.N), U_RR(conv<Mat<RR>>(U)),
          primes(primes), shift(target_shift), n(primes.length()), min_reduce_ratio(settings.reduce_ratio),
          max_level(settings.max_level), log_t(NewEnum::precomputeLogT(n)), threshold(power_ZZ(primes(n),3)),
          log_V(NewEnum::precomputeVolumes(n))
{
    // setting up the checkForEquation workspace
    this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
    this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1

    this->delayedStagesCounter.resize(this->max_level - 9, 0);     // 0 is used for total counter
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
    inv(alpha_nm1,alpha_nm1 - floor(alpha_nm1));        // this is actually alpha_n
    FloorToZZ(a_nm1,alpha_nm1);                         // this is actually a_n
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

void NewEnum::prepare(unsigned long round, const Mat<ZZ> &newBasis, const Mat<ZZ> &newU_scaled,
                      const Vec<RR> &new_target_coordinates)
{
    this->round = round;
    conv(this->U_scaled_RR,newU_scaled);
    this->tau = new_target_coordinates;

    this->current_level = 10; //this->min_level;
    this->clearL();

    ComputeGS(transpose(newBasis),this->mu, this->R_ii_squared);

    this->equations.clear();
    for(long i = 0; i < this->max_level - 10 + 1; i++)
        this->delayedStagesCounter[i] = 0;

    // setting the decreasing behavior
    this->decrease_max_distance = true;

    this->precomputeLogV();

    // start algorithm with a start parameter A
    this->A_curr = 0;
    for(long i = 1; i <= this->R_ii_squared.length(); i++)
    {
        this->A_curr += this->R_ii_squared(i);
    }

    this->A_curr *= 0.25;
    this->theoreticalMaxDistance = this->A_curr;
    this->A_curr *= this->settings.A_start_factor;
    this->heuristicMaxDistance = this->A_curr;
}