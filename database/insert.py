import database.connect 

def select_push_target(registration_table,
						results_table,
						success_data_path):
	"""
	This method is to inform job pusher which targets 
	to push, which need meet three requirements:
	1. status is job_success
	2. job files (.log and .inp) located as expected
	3. results table doesn't have this job at
	   that level of theory
	Returns a list of targets with necessary meta data
	"""
	reg_query = {"status":"job_success"}
	targets = list(registration_table.find(reg_query))

	selected_targets = []
	for target in targets:
		aug_inchi = str(target['aug_inchi'])
		spec_name = aug_inchi.replace('/', '_slash_')
		spec_path = os.path.join(success_data_path, spec_name)
		log_path = os.path.join(spec_path, 'input.log')
		inp_path = os.path.join(spec_path, 'input.inp')
		if os.path.exists(log_path) and os.path.exists(inp_path):
			level_of_theory = autoqm.utils.get_level_of_theory(inp_path)

			# query results table
			res_query = {"aug_inchi":aug_inchi, 
						"level_of_theory":level_of_theory}
			res_entries = list(results_table.find(res_query))
			if len(res_entries) == 0:
				# means no records of this target
				# in results table
				selected_targets.append(target)
	return selected_targets