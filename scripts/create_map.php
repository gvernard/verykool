<?php


$path = $argv[1]; // Absolute path to a MultiNest run of the code



# Get MAP parameters and round them at the 3rd decimal
######################################################
$json = json_decode(file_get_contents($path . 'output/smooth_output.json'),true);

$phys_maps = array();
foreach($json['json_active']['physical'] as $obj){
    $phys_maps[$obj["nam"]] = round($obj["map"],3);
}
$reg_maps = array();
foreach($json['json_active']['reg'] as $obj){
    $reg_maps[$obj["nam"]] = round($obj["map"],3);
}
$lens_maps = array();
foreach($json['json_active']['lenses']['dum0'] as $obj){
    $lens_maps[$obj["nam"]] = round($obj["map"],3);
}



# Get vkl_input.json
######################################################
$json_in = json_decode(file_get_contents($path . 'vkl_input.json'),true);

for($i=0;$i<count($json_in["physical"]["nlpars"]);$i++){
    foreach( $phys_maps as $key => $val ){
	if( $key == $json_in["physical"]["nlpars"][$i]["nam"] ){
	    $json_in["physical"]["nlpars"][$i]["val"] = $val;
	}
    }
}
for($i=0;$i<count($json_in["reg_s"]["nlpars"]);$i++){
    foreach( $reg_maps as $key => $val ){
	if( $key == $json_in["reg_s"]["nlpars"][$i]["nam"] ){
	    $json_in["reg_s"]["nlpars"][$i]["val"] = $val;
	}
    }
}
for($i=0;$i<count($json_in["lenses"]["dum0"]["nlpars"]);$i++){
    foreach( $lens_maps as $key => $val ){
	if( $key == $json_in["lenses"]["dum0"]["nlpars"][$i]["nam"] ){
	    $json_in["lenses"]["dum0"]["nlpars"][$i]["val"] = $val;
	}
    }
}
$json_in["nproc"] = 1;
$json_in["minimizer"]["type"] = "test";



# Get MAP parameters and round them at the 3rd decimal
######################################################
$dum = explode('/',trim($path,'/'));
$case = array_pop($dum);


$new_path = '/' . implode('/',$dum) . '/'. $case . "_map/";
if( !is_dir($new_path) ){
    mkdir($new_path);
}
file_put_contents($new_path . 'vkl_input.json',json_encode($json_in,JSON_PRETTY_PRINT | JSON_UNESCAPED_SLASHES));

?>
