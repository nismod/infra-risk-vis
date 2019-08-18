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

function unique(arr) {
  return Array.from(new Set(arr))
}

function insert_string(str){
	if (str.includes("bridge"))
	{
		return (str.slice(0, 10) + "to climate resilient " + str.slice(9))
	}
	else
	{
		return (str.slice(0, 13) + "climate resilient " + str.slice(13))
	}

}

export { commas, titleCase, unique, insert_string}
