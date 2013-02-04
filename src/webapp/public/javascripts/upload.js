var counter = 0;

function submit_activity() {
	$("#upload_form").collapse('hide');
	$("#upload_running").collapse('show');

	setInterval("count_up_working('upload_running','Uploading for ')", 1000); 
}

function count_up_working(id_name, status_text) {
	console.log(status_text + counter + " seconds.");
	document.getElementById(id_name).innerHTML = status_text + counter + " seconds.";
	counter = counter + 1;
}
