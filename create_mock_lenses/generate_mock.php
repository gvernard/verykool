<?php

$res_path = '';
$name = '';

if( $argc < 3 ){
    echo "\n";
    echo "==> Need command line arguments: the name of the mock observation and the name of the input file\n";
    echo "\n";
    exit;
} else {
    $name    = $argv[1];
    $options = $argv[2];

    $dum = explode('/',$options);
    if( count($dum) > 1 ){
	$tmp = array_slice($dum,0,-1);
	$res_path = implode('/',$tmp) . "/";
    } else {
	$res_path = '';
    }
}

if( file_exists($res_path . $name) ){
    echo "Output directory '" . $res_path . $name . "' exists. Overwrite? (y/n): ";
    $fh = fopen ("php://stdin","r");
    $line = fgets($fh);
    if( trim($line) != 'y' ){
	echo "ABORTING!\n";
	exit;
    } else {
	exec('rm -r ' . $res_path . $name . "/");
	exec('rm -r ' . $res_path . 'recreate_' .$name . "/");
    }
}





// STEP 1: create input files for the codes
//==========================================================================================================================================================
$recreate_dir = $res_path . 'recreate_' . $name;
if( file_exists($recreate_dir) ){
    exec('rm -r ' . $recreate_dir);
}
mkdir($recreate_dir);

mkdir($res_path . $name);
mkdir($res_path . $name . '/data');
mkdir($res_path . $name . '/base_run');

$output = $res_path . $name . '/data/';
$master_json = json_decode(preg_replace("#(/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/)|([\s\t]//.*)|(^//.*)#",'',file_get_contents($options)),true);

// create perturbations input
if( isset($master_json["perturbations"]) ){
    $pert_json = array();
    if( $master_json["perturbations"]["type"] == "grf" ){
	$pert_json["logvar"]  = $master_json["perturbations"]["logvar"];
	$pert_json["slope"]   = $master_json["perturbations"]["slope"];
	$pert_json["width"]   = $master_json["general"]["iplane"]["width"];
	$pert_json["height"]  = $master_json["general"]["iplane"]["height"];
	if( isset($master_json["perturbations"]["Nx"]) ){
	    $pert_json["Nwidth"] = $master_json["perturbations"]["Nx"];
	} else {
	    $pert_json["Nwidth"] = $master_json["general"]["iplane"]["pix_x"];
	}
	if( isset($master_json["perturbations"]["Ny"]) ){
	    $pert_json["Nheight"] = $master_json["perturbations"]["Ny"];
	} else {
	    $pert_json["Nheight"] = $master_json["general"]["iplane"]["pix_y"];
	}
	if( isset($master_json["perturbations"]["seed"]) ){
	    $pert_json["seed"] = $master_json["perturbations"]["seed"];
	}
	$pert_json["output"] = $output;
    }
    file_put_contents($recreate_dir . '/perturbations.json',json_encode($pert_json,JSON_UNESCAPED_SLASHES|JSON_PRETTY_PRINT));
}


// create fproject input
$fproject_json           = $master_json["fproject"];
$fproject_json['iplane'] = $master_json["general"]["iplane"];
if( isset($master_json["perturbations"]) ){
    $perts = array();
    $perts["type"] = "mass";
    $perts["subtype"] = "pert";
    $perts["pars"] = array();
    $perts["pars"]["filename"] = $output . "dpsi.fits";
    if( isset($master_json["perturbations"]["Nx"]) ){
	$perts["pars"]["Ni"] = $master_json["perturbations"]["Ny"];
	$perts["pars"]["Nj"] = $master_json["perturbations"]["Nx"];
    } else {
	$perts["pars"]["Ni"] = $master_json["general"]["iplane"]["pix_y"];
	$perts["pars"]["Nj"] = $master_json["general"]["iplane"]["pix_x"];	
    }
    $perts["pars"]["width"]    = $master_json["general"]["iplane"]["width"];
    $perts["pars"]["height"]   = $master_json["general"]["iplane"]["height"];
    $fproject_json["lenses"][] = $perts;
}
$fproject_json['output'] = $output;
file_put_contents($recreate_dir . '/fproject.json',json_encode($fproject_json,JSON_UNESCAPED_SLASHES|JSON_PRETTY_PRINT));


// create addmachine input
$addmachine_json            = $master_json["addmachine"];
$addmachine_json["iplane"]  = $master_json["general"]["iplane"];
$addmachine_json["noise"]   = $master_json["noise"];
$addmachine_json["imgpath"] = $output . "image_clean.fits";
$addmachine_json['output']  = $output;
file_put_contents($recreate_dir . '/addmachine.json',json_encode($addmachine_json,JSON_UNESCAPED_SLASHES|JSON_PRETTY_PRINT));


// create automask input
$automask_json               = $master_json["automask"];
$automask_json["iplane"]     = $master_json["general"]["iplane"];
$automask_json['output']     = $output;
file_put_contents($recreate_dir . '/automask.json',json_encode($automask_json,JSON_UNESCAPED_SLASHES|JSON_PRETTY_PRINT));






// STEP 2: run the codes
//==========================================================================================================================================================
$dum = explode("/",$argv[0]);
if( count($dum) > 1 ){
    $tmp = array_slice($dum,0,-1);
    $exe_path = implode("/",$tmp) . "/";
} else {
    $exe_path = "";
}

if( isset($master_json["perturbations"]) ){
    echo "Running Perturbations . . . ";
    if( $master_json["perturbations"]["type"] == "grf" ){
	exec('python '.$exe_path.'createPerturbations/generate_grf.py ' . $recreate_dir . '/perturbations.json');
    }
    exec($exe_path.'createPerturbations/kappa ' . $recreate_dir . '/perturbations.json');
    echo "done!\n";    
}

echo "Running FProject . . . ";
exec($exe_path.'FProject/fproject ' . $recreate_dir . '/fproject.json');
echo "done!\n";

echo "Running addMachine . . . ";
exec($exe_path.'addMachine/addmachine ' . $recreate_dir . '/addmachine.json');
echo "done!\n";

echo "Running autoMask . . . ";
//exec('python '.$exe_path.'autoMask/automask.py ' . $recreate_dir . '/automask.json');
exec($exe_path.'autoMask2/automask2 ' . $recreate_dir . '/automask.json');
echo "done!\n";






// STEP 3: manipulate output of the codes
//==========================================================================================================================================================
//==== Rename the image of the true source from 'vkl_source.fits' to 'source.fits'
rename($res_path . $name . '/data/vkl_source.fits',$res_path . $name . '/data/source.fits');

//==== Copy the PSF from the PSF database to the mock data directory
copy($addmachine_json['psfpath'],$res_path . $name . '/data/psf.fits');

//==== Copy the noise file from to the mock data directory
if( $master_json['noise']['noise_flag'] == 'uniform' ){
    // move file generated by addmachine
} else if( $master_json['noise']['noise_flag'] == 'map' or $master_json['noise']['noise_flag'] == 'realization' ){
    copy($master_json['noise']['noise_map'],$res_path . $name . '/data/noise.fits');
} else {
    // copy the covariance matrix file
}





// STEP 4: create vkl_input.json for the modelling code
//==========================================================================================================================================================
$input = array();

$data_constants = array(
    'physical' => array(
	'g'=>   array("per" => 0, "err" => 0, "min" => 0.00,  "max" => 0.09,   "pri" => "uni", "fix" => 0),
	'phi'=> array("per" => 1, "err" => 0, "min" => 0.0,   "max" => 180,    "pri" => "uni", "fix" => 0)
    ),
    'super' =>  array(
	'x0'=>  array("per" => 0, "err" => 0, "min" => -0.1, "max" => 0.1,   "pri" => "uni", "fix" => 0),
	'y0'=>  array("per" => 0, "err" => 0, "min" => -0.1, "max" => 0.1,   "pri" => "uni", "fix" => 0),
	'pa'=>  array("per" => 1, "err" => 0, "min" =>  0.0, "max" => 180,   "pri" => "uni", "fix" => 0)
    ),
    'sie' =>    array(
	'b'=>   array("per" => 0, "err" => 0, "min" =>  0.0, "max" => 5.0,   "pri" => "uni", "fix" => 0),
	'q'=>   array("per" => 0, "err" => 0, "min" =>  0.5, "max" => 0.999, "pri" => "uni", "fix" => 0)
    ),
    'spemd' =>  array(
	'b'=>   array("per" => 0, "err" => 0, "min" =>  0.0, "max" => 5.0,   "pri" => "uni", "fix" => 0),
	'q'=>   array("per" => 0, "err" => 0, "min" =>  0.5, "max" => 0.999, "pri" => "uni", "fix" => 0),
	'e'=>   array("per" => 0, "err" => 0, "min" =>  0.5, "max" => 3.0,   "pri" => "uni", "fix" => 0),
	's'=>   array("per" => 0, "err" => 0, "min" =>  0.0, "max" => 1.,    "pri" => "uni", "fix" => 1)
    ),
    'pert' => array(
    )
);

// paths to the data/output
$input['imgpath']  = 'data/image.fits';
if( $master_json['noise']['noise_flag'] == 'uniform' ){
    $input['noisepath']  = 'data/noise.dat';
} else if( $master_json['noise']['noise_flag'] == 'map' or $master_json['noise']['noise_flag'] == 'realization' ){
    $input['noisepath']  = 'data/noise.fits';
} else {

}
$input['maskpath'] = 'data/mask.fits';
$input['psfpath']  = 'data/psf.fits';
$input['output']   = 'output/';

// image plane
$input['iplane'] = array('pix_x' => $fproject_json['iplane']['pix_x'],
			 'pix_y' => $fproject_json['iplane']['pix_y'],
			 'siz_x' => $fproject_json['iplane']['width'],
			 'siz_y' => $fproject_json['iplane']['height']
);

// psf and noise
$input['noise_flag'] = $master_json['noise']['noise_flag'];
$input['psf'] = array('pix_x'  => $addmachine_json['psf']['pix_x'],
		      'pix_y'  => $addmachine_json['psf']['pix_y'],
		      'crop_x' => 10,
		      'crop_y' => 10
);

// source plane
$input['splane'] = array('type' => 'adaptive',
			 'mode' => 'image',
			 'spacing' => 3
);
$input['interp'] = 'bilinear';

// parameter model
$input['parameter_model'] = 'standard';

// processors
$input['nproc'] = 1;

// add perturbations if necessary
for($i=0;$i<count($fproject_json['lenses']);$i++){
    $subtype = $fproject_json['lenses'][$i]['subtype'];
    if( $subtype == 'pert' ){
	
	$input['perturbations'] = array();
	$input['perturbations']['parameter_model'] = 'perturbations_standard';
	$input['perturbations']['pix_x'] = 20;
	$input['perturbations']['pix_y'] = 20;

	$input['perturbations']['minimizer'] = array('type'    => 'test',
						     'nlive'   => 100,
						     'efr'     => 0.3,
						     'tol'     => 0.5,
						     'seed'    => 321,
						     'maxiter' => 6000
	);

	$input['perturbations']['reg_s'] = array();
	$input['perturbations']['reg_s']['type'] = 'curvature';
	$input['perturbations']['reg_s']['nlpars'] = array(
	    array('nam' => 'lambda_s',
		  'fix' => 0,
		  'per' => 0,
		  'val' => 1.0,
		  'err' => 0,
		  'min' => 0.0001,
		  'max' => 1000,
		  'pri' => array(
		      'type' => 'log10'
		  )
	    )
	);

	$input['perturbations']['reg_dpsi'] = array();
	$input['perturbations']['reg_dpsi']['type'] = 'curvature';
	$input['perturbations']['reg_dpsi']['nlpars'] = array(
	    array('nam' => 'lambda_dpsi',
		  'fix' => 0,
		  'per' => 0,
		  'val' => 1.0,
		  'err' => 0,
		  'min' => 0.0001,
		  'max' => 1000,
		  'pri' => array(
		      'type' => 'log10'
		  )
	    )
	);

	break;
    }
}

// sampling
$input['minimizer'] = array('type'    => 'test',
			    'nlive'   => 100,
			    'efr'     => 0.3,
			    'tol'     => 0.5,
			    'seed'    => 321,
			    'maxiter' => 6000
);

// regularization
$input['reg'] = array();
$input['reg']['type'] = 'curvature';
$input['reg']['nlpars'] = array(
    array('nam' => 'lambda',
	  'fix' => 0,
	  'per' => 0,
	  'val' => 1.0,
	  'err' => 0,
	  'min' => 0.0001,
	  'max' => 1000,
	  'pri' => array(
	      'type' => 'log10'
	  )
    )
);

// physical parameters (g,phi)
$input['physical'] = array();
$input['physical']['nlpars'] = array();
for($i=0;$i<count($fproject_json['physical']);$i++){
    $nam = $fproject_json['physical'][$i]['nam'];
    $input['physical']['nlpars'][] = array('nam' => $nam,
					   'fix' => $data_constants['physical'][$nam]['fix'],
					   'per' => $data_constants['physical'][$nam]['per'],
					   'val' => $fproject_json['physical'][$i]['val'],
					   'err' => $data_constants['physical'][$nam]['err'],
					   'min' => $fproject_json['physical'][$i]['min'],
					   'max' => $fproject_json['physical'][$i]['max'],
					   //'min' => $fproject_json['physical'][$i]['val'] - abs($fproject_json['physical'][$i]['val'])/10.0,
					   //'max' => $fproject_json['physical'][$i]['val'] + abs($fproject_json['physical'][$i]['val'])/10.0,
					   'pri' => array(
					       'type' => $data_constants['physical'][$nam]['pri']
					   )
    );
}
   
// lens parameters
$input['lenses'] = array();
for($i=0;$i<count($fproject_json['lenses']);$i++){
    $subtype = $fproject_json['lenses'][$i]['subtype'];

    if( $subtype != 'pert' ){
	
	$super_sub = array_merge($data_constants['super'],$data_constants[$subtype]);
	
	$nlpars = array();
	for($j=0;$j<count($fproject_json['lenses'][$i]['nlpars']);$j++){
	    $nam = $fproject_json['lenses'][$i]['nlpars'][$j]['nam'];
	    $nlpars[] = array('nam' => $nam,
			      'fix' => $super_sub[$nam]['fix'],
			      'per' => $super_sub[$nam]['per'],
			      'val' => $fproject_json['lenses'][$i]['nlpars'][$j]['val'],
			      'err' => $super_sub[$nam]['err'],
			      'min' => $fproject_json['lenses'][$i]['nlpars'][$j]['min'],
			      'max' => $fproject_json['lenses'][$i]['nlpars'][$j]['max'],
			      //'min' => $fproject_json['lenses'][$i]['nlpars'][$j]['val'] - abs($fproject_json['lenses'][$i]['nlpars'][$j]['val'])/10.0,
			      //'max' => $fproject_json['lenses'][$i]['nlpars'][$j]['val'] + abs($fproject_json['lenses'][$i]['nlpars'][$j]['val'])/10.0,
			      'pri' => array(
				  'type' => $super_sub[$nam]['pri']
			      )
	    );
	}
	
	$lens = array('type'    => 'mass',
		      'subtype' => $subtype,
		      'nlpars'  => $nlpars
	);
	
	$input['lenses']['dum' . $i] = $lens;
    }
}
file_put_contents($res_path . $name . '/base_run/vkl_input.json',json_encode($input,JSON_UNESCAPED_SLASHES | JSON_PRETTY_PRINT));






//exec('rm -r ' . $recreate_dir);
echo "Done! Mock lens in: " . $res_path . $name . "\n";
?>
