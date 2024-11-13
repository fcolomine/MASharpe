::@goto %1
:: This label performs the run of the
:: monitor, selector and variator.
:run
@shift

@cd %1_win
copy ..\%2_win\PISA_cfg PISA_cfg
start /min /b %1.exe %1_param.txt ../PISA_ 0.01
::start /min /b %1.exe %1_param.txt ..\PISA_ 0.1
@cd ../%2_win 
::start /min /b %2.exe %2_param.txt PISA_ 0.01
::@cd ../%1_win
%2.exe %2_param.txt ../PISA_ ../runs/%2_%1 %3 0.01 %4 %5
@cd ..
@goto end

:end
