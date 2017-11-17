#!/usr/bin/python

# Script for debugging features of automatic.py

import os
import abaverify as av

av.Automatic.generateRunTimePlots2(template='template_run_time_plots_compdam', 
	path_to_archived_tests=os.path.abspath('archivedTestResults'), 
	saveAs='testing_run_time_plots.html')