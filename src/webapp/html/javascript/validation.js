function validate_thresh(id_name) {
	var thresh_val=document.getElementById(id_name).value;
	if (thresh_val == null || thresh_val == "" || ! isFinite(thresh_val) || thresh_val <= 0) {
		document.getElementById(id_name).style.border='2px solid red'
		document.getElementById('thresh_help').className='error';
		document.getElementById('submit_button').disabled=true;
	} else {
		document.getElementById(id_name).style.border='1px solid #ccc'
		document.getElementById('thresh_help').className='hide';
		document.getElementById('submit_button').disabled=false;
	}
}

function validate_min_ad_size(id_name) {
	var thresh_val=document.getElementById(id_name).value;
	if (thresh_val <= 0) {
		document.getElementById(id_name).style.border='2px solid red'
		document.getElementById('min_fa_size_help').className='error';
		document.getElementById('submit_button').disabled=true;
	} else {
		document.getElementById(id_name).style.border='1px solid #ccc'
		document.getElementById('min_fa_size_help').className='hide';
		document.getElementById('submit_button').disabled=false;
	}
}

function validate_max_ad_size(id_name) {
	var thresh_val = document.getElementById(id_name).value;
	var min_fa_size = document.getElementById('min_adhesion_size').value;
	
	if (thresh_val == "") {
		document.getElementById(id_name).style.border='1px solid #ccc'
		document.getElementById('max_fa_size_help').className='hide';
		document.getElementById('submit_button').disabled=false;
	} else {
		if (thresh_val <= 0 || min_fa_size >= thresh_val) {
			document.getElementById(id_name).style.border='2px solid red'
			document.getElementById('max_fa_size_help').className='error';
			document.getElementById('submit_button').disabled=true;
		} else {
			document.getElementById(id_name).style.border='1px solid #ccc'
			document.getElementById('max_fa_size_help').className='hide';
			document.getElementById('submit_button').disabled=false;
		}
	}
}

function validate_email(id_name) {
	var emailPattern = /^[a-zA-Z0-9._-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,4}$/;
	var email_val=document.getElementById(id_name).value;
	if ((email_val != null && email_val != "") &&
		! emailPattern.test(email_val)) {
		document.getElementById(id_name).style.border='2px solid red';
		document.getElementById('email_help').className='error';
		document.getElementById('submit_button').disabled=true;
	} else {
		document.getElementById(id_name).style.border='1px solid #ccc'
		document.getElementById('email_help').className='hide';
		document.getElementById('submit_button').disabled=false;
	}
}

function validate_phone(id_name) {
	var phone_no_dash = /\d{10}/;
	var phone_dash = /\d{3}-\d{3}-\d{4}/;
	
	var phone_val = document.getElementById(id_name).value;
	if ((phone_val != null && phone_val != "") && 
		(! phone_no_dash.test(phone_val) && ! phone_dash.test(phone_val))){
		document.getElementById(id_name).style.border='2px solid red';
		document.getElementById('phone_help').className='error';
		document.getElementById('submit_button').disabled=true;
	
	} else {
		document.getElementById(id_name).style.border='1px solid #ccc'
		document.getElementById('phone_help').className='hide';
		document.getElementById('submit_button').disabled=false;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Adding Status Messages to Uploads/Processing
///////////////////////////////////////////////////////////////////////////////

var counter = 0;

//For thresholding testing submit
function threshold_submit_activity() {
	setInterval("count_up_working('working','Processing for ')", 1000);
}

//For full experimental submit
function submit_activity() {
	document.getElementById('upload_form').className='hide';
	document.getElementById('upload_working').className='';
	document.getElementById('uploading_notification').className='';
	
	setInterval("count_up_working('uploading_notification','Uploading for ')", 1000);
}

function count_up_working(id_name, status_text) {
	document.getElementById(id_name).innerHTML = status_text + counter + " seconds.";
	counter = counter + 1;
}

