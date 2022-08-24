# mpx_networkmodel_mmwr

Here, we provide code used to duplicate the results from the MMWR "Modeling the impact of sexual networks in the transmission of 
Monkeypox virus between gay, bisexual, and other men who have sex with men". 

To duplicate these results, run the file **run_network_models_mmwr.R**. This code will load in **mpx_network_object.rda**, which contains 
target statisticts for a sexual network of Gay, Bisexual, and other men who have sex with men (E.g., the proportion of individuals in the 
model that should have long term sexual partners). The code will additionally load in functions contained in **network_modules.R**. These 
functions will control the dynamics of our simulations (e.g., individuals initiating and ending sexual contact, individuals 
contracting and recovering from Monkeypox).  Finally, the code will run simulations for all scenarios presented in the MMWR. 

The code will create and save 12 files- 

**mmwr_nobehavechange_hightrans.Rdata** : Raw data for higher transmission scenario where individuals do not reduce the rate of acquiring one time sexual contacts. 
**mmwr_behavechange_hightrans.Rdata** : Raw data for higher transmission scenario where individuals reduce the rate of acquiring one time sexual contacts. 
**mmwr_nobehavechange_lowtrans.Rdata** : Raw data for lower transmission scenario where individuals do not reduce the rate of acquiring one time sexual contacts. 
**mmwr_behavechange_lowtrans.Rdata** : Raw data for lower transmission scenario where individuals reduce the rate of acquiring one time sexual contacts. 

**mmwr_nobehavechange_hightrans_noextinct_means.csv** : Means for results of higher transmission scenario where individuals do not reduce the 
rate of acquiring one time sexual contacts
**mmwr_behavechange_hightrans_noextinct_means.csv** : Means for results of higher transmission scenario where individuals reduce the rate of 
acquiring one time sexual contacts
**mmwr_nobehavechange_lowtrans_noextinct_means.csv** : Means for results of lower transmission scenario where individuals do not reduce the 
rate of acquiring one time sexual contacts
**mmwr_behavechange_lowtrans_noextinct_means.csv** : Means for results of lower transmission scenario where individuals reduce the rate of 
acquiring one time sexual contacts

**mmwr_nobehavechange_hightrans_noextinct_medians.csv** : Medians and interquartile ranges for results of higher transmission scenario where 
individuals do not reduce the rate of acquiring one time sexual contacts
**mmwr_behavechange_hightrans_noextinct_medians.csv** : Medians and interquartile ranges for results of higher transmission scenario where 
individuals reduce the rate of acquiring one time sexual contacts
**mmwr_nobehavechange_lowtrans_noextinct_medians.csv** : Medians and interquartile ranges for results of lower transmission scenario where 
individuals do not reduce the rate of acquiring one time sexual contacts
**mmwr_behavechange_lowtrans_noextinct_medians.csv** : Medians and interquartile ranges for results of lower transmission scenario where 
individuals reduce the rate of acquiring one time sexual contacts


CSV files of mean values for scenarios with no behavior change can be used to recreate figure 1 by using the **se.flow.main.mean**, 
**se.flow.pers.mean**, and **se.flow.ot.mean** variables, which are the average per day number of individuals who contracted monkeypox 
via "main" partnerships, "casual" partnerships, and "one time" partnerships respectively. 

CSV files of median and IQR values for all scenarios can be used to recreate figure 2 by using the **cuml.infs.med**, 
**cuml.infs.iqr1**, and **cuml.infs.iqr3** variables, which are the median, first quartile, and third quartile of the cumulative number of 
individuals who have contracted monkeypox over time. 
