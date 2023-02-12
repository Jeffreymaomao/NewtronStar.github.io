var longmath = document.querySelectorAll('.muti-line,.one-line');
function Switch(){
	console.log('switch');
	longmath.forEach((e)=>{
		e.classList.toggle("none");
	})
}