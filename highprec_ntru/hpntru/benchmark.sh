#!/bin/bash
P_CORE_NUM=0
/usr/local/bin/p_boost $P_CORE_NUM
SET=$(seq 1 1)
for i in $SET
do
    # taskset --cpu-list $P_CORE_NUM go test -benchmem -benchtime=1x -run ^$ -bench ^BenchmarkCompRLWENTRUP3 > p3_$i.txt
    # taskset --cpu-list $P_CORE_NUM go test -benchmem -benchtime=1x -run ^$ -bench ^BenchmarkCompRLWENTRUP4 > p4_$i.txt
    taskset --cpu-list $P_CORE_NUM go test -benchmem -benchtime=1x -run ^$ -bench ^BenchmarkNTRUTreeP4 > p4_5.txt
done
/usr/local/bin/p_cooldown $P_CORE_NUM