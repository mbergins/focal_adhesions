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
