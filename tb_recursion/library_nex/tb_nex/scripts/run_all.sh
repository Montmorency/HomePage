pushd ./exbccHrecal
./exbccHrecal.x > bccH.out
./exproc.x < fort.1 > exproc.out
python ../process_tb.py
popd

pushd ./exbccHSrecal
./exbccHSrecal.x > bccH.out
./exprocS.x < fort.1 > exproc.out
python ../process_tb.py
popd
