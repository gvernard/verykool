<?php 

$path = $argv[1];
$run  = $argv[2];

if( count($argv) == 4 ){
    $lmodel = $argv[3];
    $step = "";
    $out_path = $path . $run . "output/" . $lmodel;
} else if( count($argv) == 5 ){
    $lmodel = $argv[3];
    $step = $argv[4];
    $out_path = $path . $run . "output/" . $step . "_" . $lmodel;
} else {
    echo "Either 3 or 4 command line arguments required: path, run, lmodel, <step>\n";
    echo count($argv) . " provided, exiting!!!\n";
}

$json = json_decode(file_get_contents($out_path . "_output.json"),true);
$dum = $json["collapsed_all"];
$keys = array_keys($dum);
$pars = array_values($dum);

$fh = fopen("table_pars.txt","w");
for($i=0;$i<count($keys);$i++){
    fprintf($fh,"%25s %25.4f\n",str_replace('dum','',$keys[$i]),$pars[$i]);
}
fclose($fh);

?>
