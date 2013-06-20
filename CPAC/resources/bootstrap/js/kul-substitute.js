/*!
 * jQuery Kul Substitute
 * Copyright (c) 2009 Pablo Ziliani, Kultroom
 * Version: 2.00 (02-NOV-2009)
 * Dual licensed under the MIT and GPL licenses:
 * http://www.opensource.org/licenses/mit-license.php
 * http://www.gnu.org/licenses/gpl.html
 * Requires: jQuery v1.2.6 or later
 */
;(function($) {

/* tests:

>>> var data = {"spaces, not braces": "espacios, no llaves" }

>>> $.substitute("${ spaces, not braces }, ${ non-existent }", data)
"espacios, no llaves, "

>>> $.substitute("${ spaces, not braces }, ${ non-existent }", data, {missing: "(unknown)"})
"espacios, no llaves, (unknown)"

>>> $.substitute("${ spaces, not braces }, ${ non-existent }", data, {missing:null})
"espacios, no llaves, non-existent"

>>> $.substitute("${ spaces, not braces }, ${ non-existent }", data, {missing:null, raw: true})
"espacios, no llaves, ${ non-existent }"

>>> $.substitute("${ spaces, not braces }, ${ non-existent }", data, {missing:null})
"espacios, no llaves, non-existent"

>>> $.substitute("{{ spaces, not braces }}, {{ non-existent }}", data, {missing:null, regex: $.substitute.re_django})
"espacios, no llaves, non-existent"

*/

$.substitute = (function(template, vars, settings) {
	settings = $.extend({}, $.substitute, settings);
	return template.replace(settings.regex, function(tag, name, index, template){
		return settings.parser.apply(vars, [tag, name, index, template, settings]);
	});
});
$.extend($.substitute, {
	missing: '',
	re_django: /\{\{\s*(.+?)\s*\}\}/g,
	raw: false,
	regex: /\$\{\s*([^}]+?)\s*\}/g, // e.g.: ${ var }, ${ spaces, not braces }
	parser: (function(tag, name, index, template, settings){
		return this[name] || (settings.missing === null?(settings.raw?tag:name):settings.missing);
	})
});
})(jQuery)