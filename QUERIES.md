##Structure
####settings
Contains the settings the tests were run with.

- `ID`: An ID just for the database
- `N`: The size of the number that is going to be factorized. If N is approximately 10^14, then `N` holds the number 14
- `c`: A parameter of the program
- `dim`: A parameter of the program. Sometimes called n in the program.
- `s_max`: A parameter of the program also called "the pruning level".
- `reduce_ratio`: A parameter that is most of the time fixed to 0.8 and not widely used.
- `start_reduction`: A parameter that is also named `a` in the master thesis.
- `bkz_strong`: The block size of the strong BKZ reduction
- `bkz_slight`: The block size of the slight BKZ reduction
- `seed`: The seed the program was started with
- `cf`: An indicator if continued fractions (1), centered continued fraction (2) or none of them (0) are used.
- `scaling_type`: the scaling type that is used. This is an integer between 0 and 5.
- `file`: The file where you can see the results.

###statistics
Some statistics about the tests:

- `settings_id` Foreign key to the settings
- `runtime_total`
- `runtime_per_uniq_eqn`
- `rounds`: The number of rounds in total
- `rounds_wo_dist_red`: The number of rounds without a reduction of the distance
- `rounds_w_eqns`: The number of rounds where at least one equation was found
- `uniq_eqns`: The number of uniq equations
- `duplicate_eqns`: The number of equations found at least twice
- `total_checked_stages_rwe`: Total checked stages in rounds where equations were found (rwe - rounds with equations)
- `min_checked_stages_rwe` 
- `max_checked_stages_rwe`
- `avg_checked_stages_rwe`
- `total_checked_stages_rwoe`: Total checked stages in rounds where no equations were found (rwoe - rounds without equations)
- `min_checked_stages_rwoe`
- `max_checked_stages_rwoe`
- `avg_checked_stages_rwoe`
- `min_time_newenum_total`: Min time NewEnum took (any round)
- `max_time_newenum_total`
- `avg_time_newenum_total`
- `min_time_newenum_we`: Min time NewEnum took in rounds where equations were found
- `max_time_newenum_we`
- `avg_time_newenum_we`
- `min_time_slight_bkz`: Min time the slight BKZ reduction took
- `max_time_slight_bkz`
- `avg_time_slight_bkz`
- `max_distance_reduction`: The max distance reduction. This is a relative number (relative to the GNF of the reduced lattice basis)
- `min_distance_reduction`
- `avg_distance_reduction`
- `eqns_newenum_only`: The number of equations that were found by NewEnum. (Not in the continued fraction extension)
- `min_v_value_neo`: Minimal value of `v` used in equations from NewEnum only (neo)
- `max_v_value_neo`
- `avg_v_value_neo`
- `median_v_value_neo`
- `eqns_cf_only`: The number of equations that were found by the continued fractions extension only.
- `min_v_value_cfo` Minimal value of `v` used in equations from continued fraction only (cfo)
- `max_v_value_cfo`
- `avg_v_value_cfo`
- `median_v_value_cfo`

###rounds
- `ID`: An ID used as foreign key in the `equations` table
- `settings_id`: Foreign key to the `settings` table
- `num_primes`: the number of scaled primes
- `primes`: The scaled primes it self. This is formatted as a JSON string
- `round`: The number of the round in the corresponding test

###equations
The equations $u\equiv \pm |u-vN| mod N$ that were found during the tests.

- `rounds_id`: Foreign key to the `rounds` table
- `equation_id`
- `u`: u from the equation. Stored as a JSON string containing primes and the exponent. 
- `u_min_prime`: smallest prime in u
- `u_max_prime`: biggest prime in u
- `u_median_prime`
- `u_num_primes`: The number of different prime factors in u
- `u_max_e`: The highest exponent of a prime in u 
- `v`: v from the equation. Stored as a string, since the number can get very big.
- `v_smooth`: Indicator if v is a smooth number (1) or not (0)
- `rs`: |u-vN| from the equation. ALso stored as a JSON string
- `rs_sign`: The sign of u-vN
- `rs_min_prime`: smallest prime in |u-vN|
- `rs_max_prime` biggest prime in |u-vN|
- `rs_median_prime`
- `rs_num_primes`: The number of different prime factors in |u-vN|
- `rs_max_e`: The highest exponent of a prime in |u-vN|
- `cf`: An indicator that shows if the equation was found via continued fractions (1) or from NewEnum (0)
- `dist`: The relative distance the equation was found at.
- `level`: The level the equation was found at.

##Queries
- The number of tests and the time they took in hours grouped by `N`.

```
SELECT
   N,
   COUNT(*),
   SUM(runtime_total) / 3600
 FROM statistics
   JOIN settings ON statistics.settings_id = settings.ID
   WHERE file NOT LIKE 'special\\\\%'
 GROUP BY N ASC WITH ROLLUP;
```

- The number of distinct equations found grouped by `N`

```
SELECT COUNT(DISTINCT u, rs)
FROM equations
  JOIN rounds ON equations.rounds_id = rounds.ID
  JOIN settings ON rounds.settings_id = settings.ID
WHERE file NOT LIKE 'special\\\\%'
GROUP BY N ASC WITH ROLLUP;
```

- The number of equations found in the `c` for `N = 20` test series and the number of tests grouped by `c`

```
SELECT
  c,
  COUNT(*),
  COUNT(DISTINCT settings_id)
FROM settings
  JOIN rounds ON settings.ID = rounds.settings_id
  JOIN equations ON rounds.ID = equations.rounds_id
WHERE file LIKE '%c\\\\20%' AND equations.cf = 0
GROUP BY c;
```