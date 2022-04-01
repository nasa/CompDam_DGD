#!/usr/bin/python

# Sample script for automatically running abaverify verification tests
# Use a utility like cron or a git hook to call this script as needed

import os
import inspect
from optparse import OptionParser
import abaverify as av

import pprint
pp = pprint.PrettyPrinter(indent=4)

# Set paths relative to the location of this file
pathForThisFile = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

# Location for storing automatic test results
archive_dir = os.path.join(pathForThisFile, 'archivedTestResults')

# Initialize the automatic tester
av_auto = av.Automatic(test_directory=pathForThisFile, 
						archive_directory=archive_dir, 
						repository = {
							'name': 'CompDam_DGD',
							'remote': 'origin',
							'branch': 'dev'
						},
						tests_to_run=[],
						verbose=True,
						abaqus_cmd='abq2017')


# Run the tests
result = av_auto.run()

# PNG files in testOutpu
attach = [os.path.join(os.getcwd(), 'testOutput', x) for x in os.listdir(os.path.join(os.getcwd(), 'testOutput')) if x.endswith(".png")]

# Process the results
if result:
	av_auto.generateRunTimePlots(template='template_run_time_plots_compdam')
	
	av_auto.emailResults(recipients=["andrew.c.bergan@nasa.gov", "frank.a.leone@nasa.gov"], sender="noreply@nasa.gov", 
		template='template_email_summary', attachments=attach)


# TODO - implement below:

# # Post the results to github
# html_test_list = av_auto.generateReport('my_cool_template_file')
# html_plots = av_auto.generateRunTimePlots('a_template_for_plotting')
# av_auto.commitReportToGitHub([html_test_list, html_plots], github_info)
