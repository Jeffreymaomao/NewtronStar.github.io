const btn =  document.querySelector('#sidebar-icon');
const sidebar =  document.querySelector('#sidebar');
btn.addEventListener('click',function(){
	var state = sidebar.style.display;
	if(state=='block'){
		sidebar.style.display = 'none';
	}else{
		sidebar.style.display = 'block';
	}
	btn.classList.toggle('open');
	btn.classList.toggle('close');

})