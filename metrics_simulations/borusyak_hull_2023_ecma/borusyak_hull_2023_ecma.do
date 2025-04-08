**********************************************************************
* understanding that MA growth depends on RCT + geographic centrality
/*
This follows section 2 of their NBER long version, where I assume 
that there is a 8 x 8 = 64 grid island. Population is equal (=1)
in all grids. Distance is 1 from one grid to the next if connected, 
infinity otherwise. And cost function is tau_lkt = 2^(0.1*d_lkt), where
l is the region we care about, k is the region it's connected to, 
and t is date 1 or 2.
*/
*********************************************************************

set seed 143 

* create grid ID column where grid ID is a column+row ID, column and row IDs are x and y coordinates of their end-points
* so for ex: the bottom left square will have ID (1, 1), Top left has (1, 8), Top right has (8, 8) and bottom right has (8, 1)
set obs 64
gen x = _n
gen y = _n
gen row = ceil(x/8)
gen col = mod(x, 8)
replace col = 8 if col == 0
gen l = "(" + string(row) + ", " + string(col) + ")"
list

keep l row col 
* create a column, where each grid square can be connected in theory to any grid next to it, but not diagonal 
* so for example, (1, 1) can be connected to (1, 2) and (2, 1), but not (2, 2)
gen k1 = ""
gen k2 = ""
gen k3 = ""
gen k4 = ""

replace k1 = "(" + string(row) + ", " + string(col+1) + ")" if col < 8
replace k2 = "(" + string(row) + ", " + string(col-1) + ")" if col > 1
replace k3 = "(" + string(row+1) + ", " + string(col) + ")" if row < 8
replace k4 = "(" + string(row-1) + ", " + string(col) + ")" if row > 1

drop row col

reshape long k, i(l)
drop _j 
drop if missing(k)

* delete one of the rows from the pair where l is connected to k, and k is connected to l
* so for example, if (1, 1) is connected to (1, 2), then (1, 2) is connected to (1, 1), so delete one of them
* this is done by keeping only the first row of the pair
gen pair_id = cond(l < k, l + "_" + k, k + "_" + l)
bysort pair_id (l k): keep if _n == 1
drop pair_id


* we now have a dyadic dataset with all possible road connection, and column 1: source (l), column2: destination (k)

* create a randomization rule where 50% of the possible dyads are connected
gen treatment = (runiform() < 0.5)
su treatment

* plot the x y grid with lines connecting the chosen roads
* the road starts in the centre of a grid square, so if (1, 1) and (2, 1) are connected, the line should connect (0.5, 0.5) and (1.5, 0.5)
gen l_x = real(substr(l, 2, 1))
gen l_y = real(substr(l, 5, 1))
gen k_x = real(substr(k, 2, 1))
gen k_y = real(substr(k, 5, 1))

gen road_l_x = l_x - 0.5
gen road_l_y = l_y - 0.5
gen road_k_x = k_x - 0.5
gen road_k_y = k_y - 0.5

twoway (scatter road_l_y road_l_x if treatment == 1, msymbol(o) mcolor(black) msize(medium) mlwidth(medium) mlcolor(black)) /// 
       (scatter road_k_y road_k_x if treatment == 1, msymbol(o) mcolor(black) msize(medium) mlwidth(medium) mlcolor(black)) ///
       (pcspike road_l_y road_l_x road_k_y road_k_x if treatment == 1, lcolor(black)) ///
       , legend(off) xlabel(0(1)8) ylabel(0(1)8) 
graph export "roads.png", replace

* reintroduce duplicate pairs, so that ma can be calculated for all l's
preserve 
    gen k1 = l 
    gen l1 = k 
    drop l k 
    rename l1 l
    rename k1 k 
    gen k1_x = l_x
    gen k1_y = l_y
    gen l1_x = k_x
    gen l1_y = k_y
    gen road_k1_x = road_l_x
    gen road_k1_y = road_l_y
    gen road_l1_x = road_k_x
    gen road_l1_y = road_k_y
    drop road_l_x road_l_y road_k_x road_k_y l_x l_y k_x k_y
    rename road_k1_x road_k_x
    rename road_k1_y road_k_y
    rename road_l1_x road_l_x
    rename road_l1_y road_l_y
    rename k1_x k_x
    rename k1_y k_y
    rename l1_x l_x
    rename l1_y l_y
    tempfile reversepair
    save `reversepair'
restore 

append using `reversepair'
sort l k

* create a t = 1 distance variable = a really large number that turns tau^(-1) into the order of 10^-9
* create a t = 2 distance variable = 1 if connected, and a really large number that turns tau^(-1) into the order of 10^-9
gen d_lk1 = 1000000000
gen d_lk2 = 1 if treatment == 1
replace d_lk2 = 1000000000 if treatment == 0

* create a population variable for all l's = 1
gen P_l = 1
replace P_l = 10 if (inrange(l_x, 3, 6) & inrange(l_y, 3, 6))

* create a market access variable at t =1 and =2 as: \sum_k (\tau_lkt^(-1) * P_k), where tau_lkt = 2^(0.1*d_lkt)
gen theta = 1
gen tau_lk1 = 2^((0.1*d_lk1)*(-1*theta))
gen tau_lk2 = 2^((0.1*d_lk2)*(-1*theta))

gen ma1 = tau_lk1 * P_l
bys l: egen MA1 = total(ma1)
gen ma2 = tau_lk2 * P_l
bys l: egen MA2 = total(ma2)

* calculate MA change from t = 1 to t = 2
gen ma_change = MA2 - MA1
replace ma_change = -1*ma_change 

list l k treatment d_lk1 d_lk2 ma1 ma2 MA1 MA2 ma_change

/* heatmap */
* heatmap acts very weirdly at the edges, so adding some dummy observations
* add observations with lx and ly as (0, 1) to (0, 8); (0, 9) to (8, 9); (9, 9) to (9, 1); (9, 0) to (0, 0)
forval i = 0/9 {
    insobs 1
    replace l_x = 0 if missing(l_x)
    replace l_y = `i' if missing(l_y)
}
forval i = 1/9 {
    insobs 1
    replace l_x = `i' if missing(l_x)
    replace l_y = 9 if missing(l_y)
}
forval i = 8(-1)0 {
    insobs 1
    replace l_x = 9 if missing(l_x)
    replace l_y = `i' if missing(l_y)
}
forval i = 8(-1)1 {
    insobs 1
    replace l_x = `i' if missing(l_x)
    replace l_y = 0 if missing(l_y)
}
replace ma_change = 0 if missing(ma_change)
list l k ma_change l_x l_y treatment


* plot connected grids and market access growth
twoway (contour ma_change l_y l_x, heatmap crule(linear) clegend(off) ccuts(0(-5)-40)) ///
       (scatter l_y l_x if treatment == 1, msymbol(o) mcolor(black) msize(medium) mlwidth(medium) mlcolor(black)) /// 
       (scatter k_y k_x if treatment == 1, msymbol(o) mcolor(black) msize(medium) mlwidth(medium) mlcolor(black)) ///
       (pcspike l_y l_x k_y k_x if treatment == 1, lcolor(black)) ///
       , legend(off) xtitle("") ytitle("")  xscale(range(0.5 8.5) off) yscale(range(0.5 8.5) off) ///
       xlabel(, nogrid) ylabel(, nogrid)
graph export "ma_growth.png", replace

drop if missing(l)
keep l ma_change
duplicates drop 
save "ma_growth.dta", replace

*--------------------------------------------------------------*
* Stage 2: Simulate 1000 experiments and get expected MA change
*--------------------------------------------------------------*
clear all 

* create grid ID column where grid ID is a column+row ID, column and row IDs are x and y coordinates of their end-points
* so for ex: the bottom left square will have ID (1, 1), Top left has (1, 8), Top right has (8, 8) and bottom right has (8, 1)
set obs 64
gen x = _n
gen y = _n
gen row = ceil(x/8)
gen col = mod(x, 8)
replace col = 8 if col == 0
gen l = "(" + string(row) + ", " + string(col) + ")"
list

keep l row col 
* create a column, where each grid square can be connected in theory to any grid next to it, but not diagonal 
* so for example, (1, 1) can be connected to (1, 2) and (2, 1), but not (2, 2)
gen k1 = ""
gen k2 = ""
gen k3 = ""
gen k4 = ""

replace k1 = "(" + string(row) + ", " + string(col+1) + ")" if col < 8
replace k2 = "(" + string(row) + ", " + string(col-1) + ")" if col > 1
replace k3 = "(" + string(row+1) + ", " + string(col) + ")" if row < 8
replace k4 = "(" + string(row-1) + ", " + string(col) + ")" if row > 1

drop row col

reshape long k, i(l)
drop _j 
drop if missing(k)



* delete one of the rows from the pair where l is connected to k, and k is connected to l
* so for example, if (1, 1) is connected to (1, 2), then (1, 2) is connected to (1, 1), so delete one of them
* this is done by keeping only the first row of the pair
gen pair_id = cond(l < k, l + "_" + k, k + "_" + l)
bysort pair_id (l k): keep if _n == 1
drop pair_id


* Create a base dataset of dyads (l, k) - already created before
tempfile base_dyads
save `base_dyads', replace

* Set up a file to collect MA changes for each l in each iteration
postfile maresults str10 l float l_x l_y float ma_change float simnum using ma_sims.dta, replace

* Loop over 1000 simulations
forval i = 1/1000 {
    di "Simulation `i'"
    
    use `base_dyads', clear

    * Randomly assign treatment (connected = 1)
    gen treatment = (runiform() < 0.5)

    * Duplicate dyads to calculate MA from the perspective of l
    preserve
        gen k1 = l 
        gen l1 = k 
        drop l k 
        rename l1 l
        rename k1 k 
        tempfile reversepair
        save `reversepair'
    restore
    append using `reversepair'
    
    sort l k

    * Compute x/y coordinates (if needed)
    gen l_x = real(substr(l, 2, 1))
    gen l_y = real(substr(l, 5, 1))
    
    * create a population variable for all l's = 1
    gen P_l = 1
    replace P_l = 10 if (inrange(l_x, 3, 6) & inrange(l_y, 3, 6))
     
    * Create distance at t=1 (infinite) and t=2 (1 if connected)
    gen d_lk1 = 1000000000
    gen d_lk2 = 1 if treatment == 1
    replace d_lk2 = 1000000000 if treatment == 0

    * create a market access variable at t =1 and =2 as: \sum_k (\tau_lkt^(-1) * P_k), where tau_lkt = 2^(0.1*d_lkt)
    gen theta = 1
    gen tau_lk1 = 2^((0.1*d_lk1)*(-1*theta))
    gen tau_lk2 = 2^((0.1*d_lk2)*(-1*theta))
    
    gen ma1 = tau_lk1 * P_l
    bys l: egen MA1 = total(ma1)
    gen ma2 = tau_lk2 * P_l
    bys l: egen MA2 = total(ma2)
    gen ma_change = MA2 - MA1
    replace ma_change = -1*ma_change

    * Keep one obs per l
    keep l l_x l_y ma_change
    duplicates drop 
    count
    
    gen simnum = `i'

    quietly {
    forvalues j = 1/`=_N' {
        local ll = l[`j']
        local lx = l_x[`j']
        local ly = l_y[`j']
        local mchg = ma_change[`j']
        local snum = simnum[`j']

        post maresults ("`ll'") (`lx') (`ly') (`mchg') (`snum')
    }
    }
    
}

* Close the postfile and use the results
postclose maresults
use ma_sims.dta, clear
list if l == "(4, 4)"
* Collapse to get expected MA change per location
collapse (mean) ma_change, by(l l_x l_y)
list
/* heatmap */
* heatmap acts very weirdly at the edges, so adding some dummy observations
* add observations with lx and ly as (0, 1) to (0, 8); (0, 9) to (8, 9); (9, 9) to (9, 1); (9, 0) to (0, 0)
forval i = 0/9 {
    insobs 1
    replace l_x = 0 if missing(l_x)
    replace l_y = `i' if missing(l_y)
}
forval i = 1/9 {
    insobs 1
    replace l_x = `i' if missing(l_x)
    replace l_y = 9 if missing(l_y)
}
forval i = 8(-1)0 {
    insobs 1
    replace l_x = 9 if missing(l_x)
    replace l_y = `i' if missing(l_y)
}
forval i = 8(-1)1 {
    insobs 1
    replace l_x = `i' if missing(l_x)
    replace l_y = 0 if missing(l_y)
}
replace ma_change = 0 if missing(ma_change)

* Plot the expected market access growth as a heatmap
twoway (contour ma_change l_y l_x, heatmap crule(linear) clegend(off) ccuts(0(-2)-20)) ///
       , xtitle("") ytitle("")  xscale(range(0.5 8.5) off) yscale(range(0.5 8.5) off) ///
       xlabel(, nogrid) ylabel(, nogrid)
graph export "expected_ma_growth.png", replace

drop if missing(l)
keep l ma_change
duplicates drop
rename ma_change ma_change_exp

/* stage 3: understanding what is the omitted variable and what is the source of it? */
merge 1:1 l using ma_growth.dta, keep(match) nogen

* the heatplots are not required now, so I'll convert them to being positive numbers 
replace ma_change = -1*ma_change
replace ma_change_exp = -1*ma_change_exp

gen ma_recentered = ma_change - ma_change_exp

* generate value variable: del val = beta * ma_change + e
gen val = 1 * ma_recentered + rnormal(0, 2)

reg val ma_recentered
reg val ma_change 
reg val ma_change ma_change_exp

/* what if there was an unobserved variable? would the simulated instrument purge that? */