/*
* LaTeX IT - JavaScript to Convert Latex within an HTML page into Equations
* Copyright (C) 2009 William Bateman, 2008 Waipot Ngamsaad 

* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

var LatexIT = {
	mode : 'gif',
	imgnum : 0,
	isFirefox:false,
	init : function() {
		// We need to review the support for SVG. Latest released versions are not supporting this as they should
  //  if(document.implementation.hasFeature("http://www.w3.org/TR/SVG11/feature#BasicStructure", "1.1"))
	//	  this.mode='svg';

    // browser name
		// svg support in FireFox is not allowing two images to occur currently on the same line. 
		
		var ua = navigator.userAgent.toLowerCase(); 
		if(ua.indexOf("firefox")!=-1)
		{
      // this.isFirefox = true;
		}
  },

	pre : function(txt) {
		if ( !txt.match(/<img.*?>/i) )
		{
			//Clean code
			txt=txt.replace(/<br>/mgi,"");
			txt=txt.replace(/<br \/>/mgi,"");
			//Create img tag
		//	txt = " <img src=\"http://latex.codecogs.com/"+this.mode+".latex?"+ txt +"\" /> ";
		//	txt = " <object type=\"image/svg+xml\" width=\"20\" data=\"http://latex.codecogs.com/"+this.mode+".latex?"+ txt +"\" /> ";
		  
		  if(this.mode=='svg')  
			{
				// Best for Firefox
				if(this.isFirefox) 
  		   	txt = " <object type=\"image/svg+xml\" data=\"http://latex.codecogs.com/"+this.mode+".latex?"+ txt +"\" class=\"latex\" style=\"margin:0; padding:0; border:0\" /> ";
				else // Best for Chrome
  		   //	txt = " <object type=\"image/svg+xml\" data=\"http://latex.codecogs.com/"+this.mode+".latex?"+ txt +"\" class=\"latex\" /> ";
					txt = " <img src=\"http://latex.codecogs.com/"+this.mode+".latex?"+ txt +"\" alt=\""+ txt +"\" title=\""+ txt +"\" border=\"0\" class=\"latex\" /> ";
			}
			else 
	   	  txt = " <img src=\"http://latex.codecogs.com/"+this.mode+".latex?"+ txt +"\" alt=\""+txt+"\" border=\"0\" class=\"latex\" /> ";
		}
		return txt;
	},
	
	latex : function(txt) {
		var html, htmlinline;
		if(this.isFirefox) {
		  html=" <object type=\"image/svg+xml\" data=\"http://latex.codecogs.com/"+this.mode+".latex?$2\" class=\"latex\" /> ";
		  htmlinline=" <object type=\"image/svg+xml\" data=\"http://latex.codecogs.com/"+this.mode+".latex?\\inline $2\" class=\"latex\" /> ";
		}
		else {
		  html=" <img src=\"http://latex.codecogs.com/"+this.mode+".latex?$2\" border=\"0\" class=\"latex\" /> ";
		  htmlinline=" <img src=\"http://latex.codecogs.com/"+this.mode+".latex?\\inline $2\" border=\"0\" class=\"latex\" /> ";
		}
	
	
	  txt=txt.replace(/(^\$|[^\\]\$)(.*?[^\\])\$/gm, htmlinline);
		txt=txt.replace(/(^\\|[^\\]\\)\[(.*?[^\\])\\\]/mg," <br/>"+html+"<br/> "); 
		txt=txt.replace(/\\\$/mg,"\$"); 
		txt=txt.replace(/\\\\(\[|\])/mg,"$1");
		
		return txt;
	},	
	
	render : function(tag, latexmode) {
		var eqn = window.document.getElementsByTagName(tag);
		for (var i=0; i<eqn.length; i++) {
			if(latexmode)
			  eqn[i].innerHTML = LatexIT.latex(eqn[i].innerHTML);
			else if (eqn[i].getAttribute("lang") == "latex" || eqn[i].getAttribute("xml:lang") == "latex") 
			  eqn[i].innerHTML = LatexIT.pre(eqn[i].innerHTML);
		} 
	},

  add : function(tag, latexmode)
	{
		if(typeof(latexmode)=='undefined') latexmode=false; 
		if(window.addEventListener) 
			window.addEventListener('load', new Function('LatexIT.render("'+tag+'", '+latexmode+')'),false);
		else 
		  window.attachEvent('onload', new Function('LatexIT.render("'+tag+'", '+latexmode+')') );
	},

	scale : function(e,scale)
	{
		e.width =(e.width*scale);
		e.height=(e.height*scale);
	}
};

LatexIT.init();
LatexIT.add('*');