/*
 * jQuery Kul Substitute
 * Copyright (c) 2009 Pablo Ziliani, Kultroom
 * Version: 2.00 (02-NOV-2009)
 * Dual licensed under the MIT and GPL licenses:
 * http://www.opensource.org/licenses/mit-license.php
 * http://www.gnu.org/licenses/gpl.html
 * Requires: jQuery v1.2.6 or later
 */
;(function($){$.substitute=(function(e,f,g){g=$.extend({},$.substitute,g);return e.replace(g.regex,function(a,b,c,d){return g.parser.apply(f,[a,b,c,d,g])})});$.extend($.substitute,{missing:'',re_django:/\{\{\s*(.+?)\s*\}\}/g,raw:false,regex:/\$\{\s*([^}]+?)\s*\}/g,parser:(function(a,b,c,d,e){return this[b]||(e.missing===null?(e.raw?a:b):e.missing)})})})(jQuery);