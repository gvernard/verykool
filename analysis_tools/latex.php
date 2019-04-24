<?php
$path = trim($argv[1]);
$run  = trim($argv[2]);

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


$dum = file_get_contents($path . $run .  "vkl_input.json");
$dum = preg_replace("#(/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/)|([\s\t]//.*)|(^//.*)#", '', $dum);
$inp = json_decode($dum,true);

$dum = file_get_contents($path . $run . "output/" . $lmodel . "_initial_output.json");
$ini = json_decode($dum,true);

$dum = file_get_contents($out_path . "_output.json");
$out = json_decode($dum,true);

$use_truth = false;
if( $use_truth ){
    //set nlpars->val to the truth value
}




// 1st latex part: Title
$str = '';
$str .= '\section*{\hfil Name: '.str_replace('_',' ',$run).' \hfil}';
file_put_contents('generated_tex_files/title.tex',$str);




// 2nd latex part: True, MAP, mean and sdev parameters
$str = '';
$active_names = array_keys($out["collapsed_active"]);
$all_names    = array_keys($out["collapsed_all"]);
$fixed        = array_diff($all_names,$active_names);

$str .= '\begin{table*}';
$str .= '\caption{Values of all the parameters.}';
$str .= '\begin{center}';
$str .= '\begin{tabular}{ r r r r r}';
$str .= " " . implode(' & ',array('parameters','true value','MAP','mean$\pm 1 \sigma$')) .  " \\\\ \n";
$str .= '\hline ';

if( $lmodel == "smooth" ){
    $str .= printParams($inp["physical"]["nlpars"],$out["collapsed_active"]);
    $str .= '\hline';
    $str .= printParams($inp["reg"]["nlpars"],$out["collapsed_active"]);
    $str .= '\hline';
    foreach($inp['lenses'] as $key=>$lens){
	for($i=0;$i<count($lens["nlpars"]);$i++){
	    $lens["nlpars"][$i]["nam"] = $key . "_" . $lens["nlpars"][$i]["nam"];	}
	$str .= printParams($lens["nlpars"],$out["collapsed_active"]);
	$str .= '\hline';
    }
} else {
    $str .= printParams($inp["perturbations"]["reg_s"]["nlpars"],$out["collapsed_active"]);
    $str .= '\hline';
    $str .= printParams($inp["perturbations"]["reg_dpsi"]["nlpars"],$out["collapsed_active"]);
}

$str .= '\end{tabular}';
$str .= '\end{center}';
$str .= '\label{tab:1}';
$str .= '\end{table*}';
file_put_contents('generated_tex_files/table1.tex',$str);














// 3rd latex part: Table with the terms that contribute to the evidence calculation
$str =  '';
$str .= '\begin{table*}';
$str .= '\caption{Values of the evidence and its terms.}';
$str .= '\begin{center}';
$str .= '\begin{tabular}{ r l }';
$str .= ' term & value \\\\' . "\n";
$str .= '\hline ';

foreach($out["terms"] as $key=>$val){
    $str .= str_replace("_","-",$key) . ' & ' . round(floatval($val),4) . ' \\\\ ' . "\n";
}
/*
$str .= ' $-E_D$                                                & ' . round(floatval(- $out["terms"]["ED"]),3) . ' \\\\ ' . "\n";
$str .= ' $-E_S$                                                & ' . round(floatval(- $out["terms"]["ES"]),3) . ' \\\\ ' . "\n";
$str .= ' $-\lambda E_S$                                        & ' . round(floatval($out["pars"]["lambda"])*floatval(- $out["terms"]["ES"]),3) . ' \\\\ ' . "\n";
$str .= ' $-E_D - \lambda E_S$                                  & ' . round(floatval(- $out["terms"]["ED"]) + floatval($out["pars"]["lambda"])*floatval(- $out["terms"]["ES"]),3) . ' \\\\ ' . "\n";
$str .= ' $-\frac{1}{2} \textrm{log} (\textrm{det} A )$         & ' . round(floatval($out["terms"]["logdetA"]),3) . ' \\\\ ' . "\n";
$str .= ' $ \frac{1}{2} \textrm{log} (\textrm{det} H^TH )$      & ' . round(floatval($out["terms"]["logdetHtH"]),3) . ' \\\\ ' . "\n";
$str .= ' $ \frac{1}{2} \textrm{log} (\textrm{det} C^{-1}_D )$  & ' . round(floatval($out["terms"]["logdetC"]),3) . ' \\\\ ' . "\n";
$str .= ' $ \frac{N_{source}}{2} \textrm{log} \lambda$          & ' . round(floatval($out["terms"]["logl"]),3) . ' \\\\ ' . "\n";
$str .= ' $-\frac{N_{data}}{2} \textrm{log} 2\pi$               & ' . round(floatval($out["terms"]["log2pi"]),3) . ' \\\\ ' . "\n";
*/

$str .= '\hline';
$str .= '\end{tabular}';
$str .= '\end{center}';
$str .= '\label{tab:2}';
$str .= '\end{table*}';
file_put_contents('generated_tex_files/table2.tex',$str);






// 4th latex part: Table of other values of interest
$str =  '';
$str .= '\begin{table*}[!htb]';
$str .= '\caption{General parameters and values of interest.}';
$str .= '\begin{center}';
$str .= '\begin{tabular}{ r l l }';

$str .= '\hline';

$str .= ' $N_{data}$        & number of data pixels             & ' . $ini["Ndata"] . ' \\\\ ' . "\n";
$str .= ' $N_{mask}$        & number of data pixels in the mask & ' . $ini["Nmask"] . ' \\\\ ' . "\n";
$str .= ' $N_{source}$      & number of source pixels           & ' . $ini["Nsource"] . ' \\\\ ' . "\n";
$str .= ' $N_{source mask}$ & number of source pixels in mask   & ' . $ini["Nsource_mask"] . ' \\\\ ' . "\n";
if( $inp['noise_flag'] == 'uniform' ){
    $cov = trim(file_get_contents($path . "data/noise.dat"));
    $str .= ' $C_D$        & covariance of the data            & uniform(' . $cov . ') \\\\ ' . "\n";
} else if( $inp['noise_flag'] == 'map' ){
    $str .= ' $C_D$        & covariance of the data            & map \\\\ ' . "\n";
}
$str .= ' $P_{data}$   & size of data pixels               & ' . sprintf("%10.4f",$ini["Psize"]) . ' \\\\ ' . "\n";
$str .= '\hline ';

if( $lmodel == "pert" ){
    $str .= ' $N_{pert}$       & number of $\delta\psi$ pixels             & ' . $ini["Npert"] . ' \\\\ ' . "\n";
    $str .= ' $N_{pert-mask}$  & number of $\delta\psi$ pixels in the mask & ' . $ini["Npert_mask"] . ' \\\\ ' . "\n";
    $str .= ' $P_{pert}$       & size of $\delta\psi$ pixels               & ' . sprintf("%10.4f",$ini["Ppert_size"]) . ' \\\\ ' . "\n";
}

if( $lmodel == "smooth" and $inp["minimizer"]["type"] == "multinest"){
    $dum = file_get_contents($out_path . '_minimizer_output.json');
    $min = json_decode($dum,true);
    $str .= ' $N_{iter}$   & number of MultiNest iterations    & ' . $min['total_samples'] . ' \\\\ ' . "\n";
    $str .= ' $N_{repl}$   & number of live point replacements & ' . $min['replacements'] . ' \\\\ ' . "\n";
}
if( $lmodel == "pert" and $inp["perturbations"]["minimizer"]["type"] == "multinest"){
    $dum = file_get_contents($out_path . '_minimizer_output.json');
    $min = json_decode($dum,true);
    $str .= ' $N_{iter}$   & number of MultiNest iterations    & ' . $min['total_samples'] . ' \\\\ ' . "\n";
    $str .= ' $N_{repl}$   & number of live point replacements & ' . $min['replacements'] . ' \\\\ ' . "\n";
}
if( $lmodel == "pert" and $inp["minimizer"]["type"] == "iterator"){
    $dum = file_get_contents($out_path . '_minimizer_output.json');
    $min = json_decode($dum,true);
    $str .= ' $N_{iter}$   & number of iterations    & ' . $min['itarations'] . ' \\\\ ' . "\n";
}

$str .= '\hline';
$str .= '\end{tabular}';
$str .= '\end{center}';
$str .= '\label{tab:3}';
$str .= '\end{table*}';
file_put_contents('generated_tex_files/table3.tex',$str);






function printParams($nlpars,$collapsed_active){
    $str = "";
    foreach($nlpars as $param){
	if( $param['fix'] == 0 ){
	    $name = $param["nam"];
	    $truth = $param["val"];
	    $map          =        round($collapsed_active[$name]["map"],4);
	    $mean_limits  =  "$" . round($collapsed_active[$name]["mean"],4);
	    $mean_limits .= "_{" . round($collapsed_active[$name]["mean"] - $collapsed_active[$name]["s1_low"],4) . "}";
	    $mean_limits .= "^{" . round($collapsed_active[$name]["s1_high"] - $collapsed_active[$name]["mean"],4) . "}$";
	    $sigma_range  = round($collapsed_active[$name]["s1_low"],4) . " - " . round($collapsed_active[$name]["s1_high"],4);
	} else {
	    $name = $param["nam"];
	    $truth = $param["val"];
	    $map = "-";
	    $mean_limits = "-";
	    $sigma_range = "-";
	}
	$str .= " " . implode(' & ',array(str_replace("_","-",$name),$truth,$map,$mean_limits,$sigma_range)) . " \\\\ \n";
    }
    return $str;
}
	
?>
