#!/usr/bin/perl

for ($processors = 1; $processors <= 1; $processors++) 
{
  for ($elements=20; $elements <= 20; $elements+=2)
  {
    $cmd = "mpiexec -n $processors bin/x86_64-linux/AnalyticLaplaceExample-debug  $elements $elements 0 $processors 1";
    $starttime=time();
    system $cmd;
    $usertime=time()-$starttime;
    print "Runs on $processors processors and $elements elements takes $usertime seconds";
  }
}
