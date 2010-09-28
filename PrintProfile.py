#!/usr/bin/env python

import pstats

stats = pstats.Stats('BinaryProfile')

print 'Total Calls'
stats.sort_stats('calls')
stats.print_stats(20)

print 'Cumulative Time'
stats.sort_stats('cumulative')
stats.print_stats(20)

print 'Total Time'
stats.sort_stats('time')
stats.print_stats(20)

