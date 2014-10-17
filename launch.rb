#!/usr/bin/env ruby

$alpha = 0.8
$alpha_inc = 0.1

while $alpha < 0.95 do
  job = fork do
    $alpha_end = $alpha + 2 * $alpha_inc
    $tmpfile = "optim_" + $alpha.to_s +".log"
    exec "nice -n 20 ./optimization #$alpha #$alpha_end 1 > #$tmpfile"
  end
  Process.detach(job)
  $alpha += $alpha_inc
end

