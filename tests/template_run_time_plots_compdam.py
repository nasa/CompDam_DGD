# Specify either 'commit' (to label the x-axis with commit sha) or 'date' (to label the x-axis with the date the test was run)
x_axis_qty = 'date'

# Plot dims in pixels
plot_height = '300'
plot_width = ''

# For grouping tests together on same axes
test_group_prefixes = ['test_C3D8R_failureEnvelope_sig11sig22', 
					'test_S4R_failureEnvelope_sig11sig22', 
					'test_C3D8R_failureEnvelope_sig12sig22',
					'test_S4R_failureEnvelope_sig12sig22',
					'test_C3D8R_failureEnvelope_sig12sig23',
					'test_C3D8R_mixedModeMatrix',
					'test_C3D8R_fiberCompression_DGD_wkbToTotal',
					'test_C3D8R_twoElement_fiberCompression_DGD']

# Grouping charts into subsections
chart_groups = dict()
chart_groups['failure-envelopes-C3D8R'] = {
	'name_pretty': 'Failure Envelopes, C3D8R', 
	'charts': ['test_C3D8R_failureEnvelope_sig11sig22', 'test_C3D8R_failureEnvelope_sig12sig22', 'test_C3D8R_failureEnvelope_sig12sig23']
}
chart_groups['failure-envelopes-S4R'] = {
	'name_pretty': 'Failure Envelopes, S4R',
	'charts': ['test_S4R_failureEnvelope_sig11sig22', 'test_S4R_failureEnvelope_sig12sig22']
}
chart_groups['single-element-C3D8R'] = {
	'name_pretty': 'Single Element, C3D8R',
	'charts': ['test_C3D8R_elastic_fiberTension', 'test_C3D8R_elastic_matrixTension', 'test_C3D8R_elementSize', 'test_C3D8R_fiberCompression_CDM', 'test_C3D8R_fiberLoadReversal', 'test_C3D8R_fiberTension', 'test_C3D8R_matrixCompression', 'test_C3D8R_matrixTension', 'test_C3D8R_nonlinearShear12', 'test_C3D8R_schapery12', 'test_C3D8R_simpleShear12', 'test_C3D8R_simpleShear12friction']
}
chart_groups['Fiber-Compression-DGD'] = {
	'name_pretty': 'Fiber Compression, DGD',
	'charts': ['test_C3D8R_fiberCompression_DGD', 'test_C3D8R_twoElement_fiberCompression_DGD']
}
chart_groups['Single_Element_S4R'] = {
	'name_pretty': 'Single Element, S4R',
	'charts': ['test_S4R_elementSize', 'test_S4R_fiberCompression_CDM', 'test_S4R_fiberLoadReversal', 'test_S4R_fiberTension', 'test_S4R_matrixCompression', 'test_S4R_matrixTension', 'test_S4R_nonlinearShear12', 'test_S4R_simpleShear12', 'test_S4R_simpleShear12friction']
}



# Subsection heading
subsection = """
<div id="{section_name_dashes}" class="">
  <h3>{section_name}</h3>
  {plots}
</div>
"""

# Subsection toc wrapper
subsection_toc_wrapper = """
<li>
  <a href="#{section_name_dashes}" class="">{section_name}</a>
  <ul class="nav subsection_heading">
  {toc_entries}
  </ul>
</li>
"""

# Formatting for each plt
plot = """
<div id="{plot_title}" class="section scrollspy">
{plot}
</div>
<br><br>
"""

# Table of contents
toc = """
<li>
  <a href="#{plot_title}" class="">{plot_title}</a>
</li>
"""

# Overall page formatting
body = ""
with open('template_run_time_plots_body.html') as f:
	body = f.read()