var i = 1;
function w3_open() {
  document.getElementById("mySidebar").style.display = "block";
}
 
function w3_close(p) {
  document.getElementById("mySidebar").style.display = "none";
  i += p;
}

function goBack() {
    window.history.go(-i);
}