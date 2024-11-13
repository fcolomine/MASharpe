for /L %%j IN (1,1,2500) DO (
nsga2.exe nsga2_param.txt ..\PISA_ 0.01
)