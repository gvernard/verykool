<?php 

$path = $argv[1];
$case = $argv[2];

if( count($argv) > 3 ){
    $step = $argv[3] . "_";
} else {
    $step = "";
}

$json = json_decode(file_get_contents($path . $case . "output/".$step."vkl_output.json"),true);
$dum = $json["full_pars"];
$keys = array_keys($dum);
$pars = array_values($dum);

$fh = fopen("table_pars.txt","w");
for($i=0;$i<count($keys);$i++){
    fprintf($fh,"%25s %25.4f\n",str_replace('dum','',$keys[$i]),$pars[$i]);
}
fclose($fh);

?>
