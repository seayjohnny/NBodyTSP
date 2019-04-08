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
  if(loc != null) {
    var element_to_scroll_to = document.getElementById(window.location.hash);
    element_to_scroll_to.scrollIntoView();
  }
}

function showTab(body, tab) {
  var elements = document.getElementsByClassName("tab-body");
  for (var i = 0; i < elements.length; i++) {
    elements[i].hidden = true;
  }
  elements = document.getElementsByClassName("tab-button");
  for (var i = 0; i < elements.length; i++) {
    elements[i].classList.toggle("selected", false);
  }

  document.getElementById(body).hidden = false;
  tab.classList.toggle("selected", true);
}