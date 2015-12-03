
#include "NewEnum.h"

void NewEnumStage::get_u(Vec<RR> &u) const
{
    conv(u,this->u);
}

bool StageStorage::storeStage(const RR& y_t, const RR& c_t, const RR& c_tp1, const Vec<RR>& u, long t, long level)
{
    NewEnumStage* stage = this->getStage();
    stage->set(y_t,c_tp1,u,t);

    double alpha_2;
    conv(alpha_2, c_t / this->maxDistance);
    long alpha_2_indicator = this->alpha_2_indicator(alpha_2);
    long t_indicator = this->t_indicator(stage->t);
    this->storage[t_indicator][alpha_2_indicator][level - this->min_level - 1].push_back(stage);
    if(this->alpha_2_min[alpha_2_indicator][t_indicator][level - this->min_level - 1] > alpha_2)
        this->alpha_2_min[alpha_2_indicator][t_indicator][level - this->min_level - 1] = alpha_2;
    this->stageCounterTotal++;
    this->stageCounterByLevel[level - this->min_level - 1]++;

    // statistics
    this->delayedStages[alpha_2_indicator][t_indicator][level - this->min_level - 1]++;

    return this->stageCounterTotal >= this->stageCounterThreshold;      // about 6 GB of RAM depending on the dimension
}

bool StageStorage::getNext(NewEnumStage* &stage)
{
    if(this->stageCounterByLevel[this->currentLevel - this->min_level - 1] <= 0 || this->stageCounterTotal >= this->stageCounterThreshold)
        return false;

    for(unsigned long t_indicator = 0; t_indicator < 3; ++t_indicator)
    {
        for(unsigned long alpha_2_indicator = 0; alpha_2_indicator < 3; ++alpha_2_indicator)
        {
            if(!this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].empty())
            {
                stage = this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].front();
                this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].pop_front();
                this->stageCounterByLevel[this->currentLevel - this->min_level-1]--;
                this->stageCounterTotal--;
                return true;
            }
        }
    }

    // It never should come this far, but who knows...
    return false;
}

void StageStorage::incrementCurrentLevel()
{
    this->currentLevel++;

    for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
        for(long t_indicator = 0; t_indicator < 3; t_indicator++)
        {
            this->maxDelayedAndPerformedStages[alpha_2_indicator][t_indicator][this->currentLevel - this->min_level - 1] = this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].size();
            this->totalDelayedAndPerformedStages += this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].size();
        }
}

void StageStorage::updateMaxDistance(const RR &distance)
{
    double alpha_1;
    conv(alpha_1, distance/this->maxDistance);

    if(alpha_1 < this->alpha_1_threshold && this->stageCounterTotal > 100)
        this->recalculateLevels(alpha_1);

    for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
        for(long t_indicator = 0; t_indicator < 3; t_indicator++)
            for(long s = std::max(0L,this->currentLevel-this->min_level-1); s < this->pruningLevel-this->min_level;s++)
                this->alpha_2_min[alpha_2_indicator][t_indicator][s] /= alpha_1;

    this->maxDistance = distance;
}

void StageStorage::recalculateLevels(const double &alpha_1)
{
    double new_alpha_2;
    long new_alpha_2_indicator;
    long level_change;
    long s_max = this->pruningLevel - this->min_level - 1;
    long s_min = this->currentLevel-this->min_level - 1;

    /**
     * skip t_indicator = 0, since the stages here have t between 4 and 20.
     * So in almost all cases the will be no change.
     */
    for(int t_indicator = 2; t_indicator >= 1; t_indicator--)
    {
        for(int alpha_2_indicator = 2; alpha_2_indicator >= 0; alpha_2_indicator--)
        {
            for(long s = s_max; s > s_min; s--)
            {
                if (this->storage[t_indicator][alpha_2_indicator][s].size() < 10)
                    continue;

                level_change = this->levelChange(alpha_1,
                                                 this->alpha_2_min[alpha_2_indicator][t_indicator][s],
                                                 t_indicator);
                if (level_change <= 0)
                    continue;

                new_alpha_2 = this->alpha_2_min[alpha_2_indicator][t_indicator][s] / alpha_1;
                new_alpha_2_indicator = this->alpha_2_indicator(new_alpha_2);

                if(s+level_change > this->pruningLevel-11)
                {   // return stages
                    this->stageCounterByLevel[s] -= this->storage[t_indicator][alpha_2_indicator][s].size();
                    this->stageCounterTotal -= this->storage[t_indicator][alpha_2_indicator][s].size();
                    this->returnStages(this->storage[t_indicator][alpha_2_indicator][s]);
                    this->alpha_2_min[t_indicator][alpha_2_indicator][s] = 2.0;
                }
                else
                {   // move stages
                    this->stageCounterByLevel[s] -= this->storage[t_indicator][alpha_2_indicator][s].size();
                    this->stageCounterByLevel[s + level_change] += this->storage[t_indicator][alpha_2_indicator][s].size();
                    this->storage[t_indicator][new_alpha_2_indicator][s + level_change].splice(
                            this->storage[t_indicator][new_alpha_2_indicator][s + level_change].end(),
                            this->storage[t_indicator][alpha_2_indicator][s]);
                    this->alpha_2_min[t_indicator][alpha_2_indicator][s + level_change] = std::min(this->alpha_2_min[t_indicator][alpha_2_indicator][s + level_change],
                                                                                                   this->alpha_2_min[t_indicator][alpha_2_indicator][s]);
                    this->alpha_2_min[t_indicator][alpha_2_indicator][s] = 2.0;
                }
            }
        }
    }
}

unsigned long StageStorage::levelChange(const double &alpha_1, const double &alpha_2, long t_indicator)
{
    return t_indicator == 0
           ? 0
           : (alpha_1 <= alpha_2
              ? (unsigned long) this->pruningLevel - this->min_level - 1        // 2*log(2)= 1.386294361119890..., 39=40-1
              : (unsigned long) ceil((t_indicator == 2 ? 39.0 : 17.0)  * 1.38629436111989 / log((1.0-alpha_2) / (alpha_1 - alpha_2))) - 1);
}

Vec<double> NewEnum::precomputeLogT(long n)
{
    Vec<double> log_t;
    log_t.SetLength(n+1);
    for(long t = 1; t <= n+1; t++)
        log_t(t) = log(t);

    return log_t;
}

void NewEnum::closest_RR(RR &out, const RR &x)
{
    sub(out,x,0.5);
    ceil(out, out);
}

void NewEnum::next(RR &out, const RR &u, const RR &y)
{
    RR temp;
    this->closest_RR(temp, y);

    float side1 = temp >= y ? 1.0f : -1.0f;
    float side2 = u >= y ? 1.0f : -1.0f;
    bool plusOneSwitch = (side1 == side2 || temp == u);

    // out = 2 * closest_y - u;
    mul(temp, temp, 2.0);
    sub(out, temp, u);

    if(plusOneSwitch)
        out -= side2;
}

void NewEnum::run(unsigned long round, const Mat<ZZ> &newBasis_transposed, const Mat<ZZ> &new_U_scaled,
                  const Vec<RR> &new_target_coordinates)
{
    this->prepare(round, newBasis_transposed, new_U_scaled, new_target_coordinates);

    cout << "Performing: ";
    // Reset list of delayed stages
    this->current_level = this->min_level;

    Vec<RR> u;                              // coordinates of a close vector
    u.SetLength(this->n);
    this->closest_RR(u(this->n),this->tau(this->n));

    NewEnumStage* stage = new NewEnumStage(this->tau(this->n), conv<RR>(0), u, this->n);

    this->perform(stage);
    this->L.returnStage(stage);

    for(long l = 0; l <= this->max_level - this->min_level - 1; l++)
    {
        this->L.incrementCurrentLevel();
        this->current_level++;

        // perform delayed stages
        while(this->L.getNext(stage))
        {
            this->perform(stage);
            this->L.returnStage(stage);
        }

        cout << ".";
    }

    cout << endl;
}

void NewEnum::perform(NewEnumStage* current_stage)
{
    this->decrease_max_distance = true;

    long max_t = current_stage->t;
    long t = current_stage->t;          // projection to the n-t-1 last coordinates
    Vec<RR> c,                          // Length in the projection
            u,                          // coordinates of a close vector
            y;                          // uncovered parts of the target vector
    c.SetLength(t+1);
    c(t+1) = current_stage->c_tp1;
    current_stage->get_u(u);
    y.SetLength(t);
    y(t) = current_stage->y_t;

    long level;

    bool success;

    RR temp;

    // perform stages with s = current_level
    while(t <= max_t)
    {
        // c(t) = c(t+1) + power(abs(u(t) - y(t)), 2) * this->R_ii_squared(t);
        sub(c(t),u(t),y(t));
        abs(c(t),c(t));
        power(c(t),c(t),2);
        c(t) *= this->R_ii_squared(t);
        c(t) += c(t+1);

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
                    this->L.updateMaxDistance(this->A_curr);
                }
            }

            this->decrease_max_distance = this->decrease_max_distance && !success;

            // 2.1
            goto cleanup;
        }

        // compute the probability beta
        if(t < 4 || (level = calculateLevel(t,c(t))) <= this->current_level)
        {
            t--;

            clear(y(t));        // y(t) = 0
            for(long i = t + 1; i <= this->n; i++)    // 1/r_tt * sum_{i=t+1}^n (\tau_i - u_i) r_ti
            {
                // y(t) += (this->tau(i) - u(i)) * this->mu(i, t);
                sub(temp,this->tau(i),u(i));
                mul(temp,temp,this->mu(i, t));
                y(t) += temp;
            }
            y(t) += this->tau(t);

            this->closest_RR(u(t),y(t));
            continue;
        }

        // the program comes only this far, if level > current_s holds
        if(level <= this->max_level)
        {
            if(this->L.storeStage(y(t),c(t),c(t+1),u,t,level))      // store the stage for later
                return;
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

    this->log_V_minus_log_R_prod_minus_log_t = this->log_V;

    for(long i = 1; i <= this->n; i++)
    {
        this->log_V_minus_log_R_prod_minus_log_t(i) -= log_R_diag_prod(i);
        this->log_V_minus_log_R_prod_minus_log_t(i) -= this->log_t(i + 1);
    }
}

bool NewEnum::checkForEquation(const Vec<RR> &input, const RR &c_1)
{
    this->stagesCheckedForEquations++;
    mul(this->temp_vec, this->U_scaled_RR, input);
    this->temp_vec += this->shift;
    mul(this->close_vec, this->U_RR, this->temp_vec);

    conv(this->close_vec_long,this->close_vec);

    this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
    this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1

    NTL::set(this->u);       // u = 1

    ZZ temp_ZZ,temp_ZZ_2;
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

    conv(this->v,this->u); conv(temp_RR,this->N);
    this->v /= temp_RR;

    this->closest_RR(temp_RR,this->v);

    sub(this->d,this->v,temp_RR);

    NTL::abs(this->alpha_nm1,this->d);

    this->closest_RR(temp_RR,this->v);
    conv(this->vN,temp_RR);
    this->vN *= this->N;

    NTL::clear(this->h_n); NTL::set(this->h_nm1); NTL::clear(this->h_nm2);
    NTL::set(this->k_n); NTL::clear(this->k_nm1); NTL::set(this->k_nm2);
    NTL::clear(this->a_nm1);

    long sign = NTL::sign(this->d),
         equation_counter = 0;
    bool cf_equation = false;

    do
    {
        this->nextContinuedFraction(this->h_n, this->k_n, this->h_nm1, this->k_nm1, this->h_nm2, this->k_nm2, this->a_nm1, this->alpha_nm1);

        this->equation = this->raw_equation;
//        this->left_side = this->u * this->k_n;
        mul(this->left_side,this->u,this->k_n);
//        abs(this->right_side,this->left_side - this->vN * this->k_n - sign * this->h_n * this->N);
        mul(temp_ZZ,this->vN,this->k_n);
        sub(temp_ZZ,this->left_side,temp_ZZ);
        mul(temp_ZZ_2,sign,this->h_n);
        mul(temp_ZZ_2,temp_ZZ_2,this->N);
        sub(this->right_side,temp_ZZ,temp_ZZ_2);
        abs(this->right_side,this->right_side);

        if(this->isSmooth(this->equation, this->k_n, this->left_side, this->right_side))
        {
            // this->ride_side_factor = conv<ZZ>(this->closest_RR(conv<RR>(left_side) / this->N_RR));
            equation_counter++;
            conv(temp_RR,left_side);
            div(temp_RR,temp_RR,this->N_RR);
            this->closest_RR(temp_RR,temp_RR);
            conv(this->ride_side_factor,temp_RR);
            div(temp_RR,c_1,this->theoreticalMaxDistance);
            this->equations.emplace_back(this->equation, this->ride_side_factor, this->current_level, conv<double>(temp_RR), this->round, this->timer.step(),cf_equation);
        }

        cf_equation = true;
    }
    while(this->k_n < this->threshold && this->settings.useContinuedFractions);

    return equation_counter > 0;
}

NewEnum::NewEnum(const FactoringSettings &settings, Timer &timer, const Vec<long> &primes,
                 const Mat<ZZ> &U, const Vec<RR> &target_shift)
        : settings(settings), timer(timer), N(settings.N), N_RR(conv<RR>(settings.N)), U_RR(conv<Mat<RR>>(U)),
          primes(primes), shift(target_shift), n(primes.length()), min_reduce_ratio(settings.reduce_ratio),
          max_level(settings.max_level), log_t(NewEnum::precomputeLogT(n)), threshold(power_ZZ(primes(n),3)),
          log_V(NewEnum::precomputeVolumes(n)), L(StageStorage((unsigned long) settings.max_level))
{
    // setting up the checkForEquation workspace
    this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
    this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1
}

list<Equation> NewEnum::getEquations()
{
    return this->equations;
}

void NewEnum::nextContinuedFraction(ZZ &h_n, ZZ &k_n, ZZ &h_nm1, ZZ &k_nm1, ZZ &h_nm2, ZZ &k_nm2, ZZ &a_nm1, RR &alpha_nm1)
{
//    h_n = a_nm1 * h_nm1 + h_nm2;
    mul(h_n,a_nm1,h_nm1);
    add(h_n,h_n,h_nm2);

//    k_n = a_nm1 * k_nm1 + k_nm2;
    mul(k_n,a_nm1,k_nm1);
    add(k_n,k_n,k_nm2);

    swap(h_nm2,h_nm1);
    h_nm1 = h_n;

    swap(k_nm2,k_nm1);
    k_nm1 = k_n;

    sub(alpha_nm1,alpha_nm1,floor(alpha_nm1));
    inv(alpha_nm1,alpha_nm1);                           // this is actually alpha_n
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
        equation(this->n + 1) = 0;

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

void NewEnum::prepare(unsigned long round, const Mat<ZZ> &newBasis_transposed, const Mat<ZZ> &newU_scaled,
                      const Vec<RR> &new_target_coordinates)
{
    this->round = round;
    conv(this->U_scaled_RR,newU_scaled);
    this->tau = new_target_coordinates;

    this->current_level = this->min_level;
    this->stagesCheckedForEquations = 0;

    ComputeGS(newBasis_transposed, this->mu, this->R_ii_squared);

    this->equations.clear();

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

    this->L.resetStorage(this->A_curr);
}

void NewEnum::getDistances(RR &theoretical,RR &heuristic,RR &reduced)
{
    theoretical = this->theoreticalMaxDistance;
    heuristic = this->heuristicMaxDistance;
    reduced = this->A_curr;
}