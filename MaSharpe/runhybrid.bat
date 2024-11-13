
:: The following variables contain the
:: the names of the selectors and variators
:: to be used. For each combination, there
:: will be runs of selector-variator pairs. 
:: Note that all directories and
:: parameter files need to be named
:: consistently. If a variator is used with 
:: different parameter files, it should be 
:: copied into a new directory with a
:: different name
:: spea2 ibea
SET selectors= ibea 
SET variators= masharpe

:: set modelo:
:: Use MARKOWITZ O SHARPE

SET modelo= MARKOWITZ

:: ============================================
:: DONT CHANGE BELOW THIS LINE
:: ============================================
:: This line calls the batch file compute.bat to
:: run the selector, variator and monitor.
:: The resulting files are put into the directory
:: ./runs. THIS DIRECTORY SHOULD BE EMPTY INITIALLY.
:: OTHERWISE, THE STATISTICAL TEST WILL BE WRONG.
:: for /L %%j IN (50,50,50) (number of generations:You can configure this cycle with the number of generations you want: 
:: (initial generation, increment, final generation)

for /L %%j IN (50,50,50) DO (

@FOR %%m IN (%modelo%) DO (
for /L %%k IN (1,1,1) DO (
@FOR %%s IN (%selectors%) DO (@FOR %%v IN (%variators%) DO @call compute run %%s %%v %%m %%j %%k)
)
)
)
pause
