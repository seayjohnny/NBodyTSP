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

function onFrameLoad() {
  var loc = window.location.hash;
  var element_to_scroll_to = document.getElementById(window.location.hash);
  element_to_scroll_to.scrollIntoView();
}