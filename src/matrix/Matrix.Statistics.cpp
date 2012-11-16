static void condensedMatrix_calculate_statistics(CondensedMatrix* self, PyObject *args){
	self->statistics_already_calculated = false;
}

static PyObject* condensedMatrix_get_mean(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->mean);
}

static PyObject* condensedMatrix_get_variance(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->variance);
}

static PyObject* condensedMatrix_get_skewness(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->skewness);
}

static PyObject* condensedMatrix_get_kurtosis(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->kurtosis);
}

static PyObject* condensedMatrix_get_max(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("f", self->statisticsCalculator->max);
}

static PyObject* condensedMatrix_get_min(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("f", self->statisticsCalculator->min);
}
