/**
 * Common helper functions
 *
 */

function commas(x) {
  return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

function titleCase(str) {
 var splitStr = str.toLowerCase().split(' ');
 for (var k = 0; k < splitStr.length; k++) {
     splitStr[k] = splitStr[k].charAt(0).toUpperCase() + splitStr[k].substring(1);
 }
 return splitStr.join(' ');
}

export { commas, titleCase }
