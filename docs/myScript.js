var i = 1;
function w3_open() {
  document.getElementById("mySidebar").style.display = "block";
}
 
function w3_close(p) {
  document.getElementById("mySidebar").style.display = "none";
  if(p==undefined){ p = 1;}
  i += p;
}

function goBack() {
    window.history.go(-i);
}
