<?php 

$html = '<html>';
$html .= '<head>';
$html .= '<meta charset="UTF-8">';
$html .= '<link rel="stylesheet" href="css/index.css">';
$html .= '<script src="https://code.jquery.com/jquery-1.12.4.js"></script>';
$html .= '<script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>';
$html .= '<script type="text/javascript" src="js/index.js"></script>';
$html .= '</head>';

$html .= '<body>';

$html .= '<h2>Creating lenses</h2>';


$html .= '<form id="myform">';
$html .= '<div id="tabs" class="container">';

$html .= '<ul class="tabs">';
$html .= '<li class="tab-link" data-tab="source_div">Source</li>';
$html .= '<li class="tab-link current" data-tab="lens_div">Lens</li>';
$html .= '<li class="tab-link" data-tab="observer_div">Observer</li>';
$html .= '</ul>';


//////////////////////////////
// source options
$html .= '<div id="source_div" class="tab-content">';
$json = json_decode(file_get_contents('form_json/options_source.json'),true);

$html .= '<div class="slot">';
$html .= '<p class="title">Source</p>';


$str = '<select class="myselect" id="source">';
$str .= '<option disabled selected="selected">Select option</option>';
foreach($json as $key => $val){
    $str .= '<option value="'.$key.'">'.$val['label'].'</option>';
}
$str .= '</select>';

// Analytic
$str .= '<div class="myselect_option analytic">';
$str .= '<p class="legend">' . $json['analytic']['description'] . '</p>';


$str .= '<hr class="myhr">';
$str .= '<div id="anal-source-1" class="anal-source">';
$str .= '<h4>Analytic Source 1</h4>';
$str .= '<button type="button" class="remove_source">Remove</button>';
$str .= createSelect($json['analytic']['select_list'],'analytic');
$str .= '<hr class="myhr">';
$str .= '</div>';
$str .= '<button type="button" id="add_source">Add source</button>';
$str .= '</div>';

// Fromfits
$str .= '<div class="myselect_option fromfits">';
$str .= '<p class="legend">' . $json['fromfits']['description'] . '</p>';
$str .= createInputs($json['fromfits']['input_list'],$json['fromfits']['type'],'fromfits');
$str .= '</div>';

// Delaunay
$str .= '<div class="myselect_option delaunay">';
$str .= '<p class="legend">' . $json['delaunay']['description'] . '</p>';
$str .= createInputs($json['delaunay']['input_list'],$json['delaunay']['type'],'delaunay');
$str .= '</div>';


$html .= $str;
#$html .= createSelect($json['sources']['select_list']);


$html .= '</div>';

$html .= '</div>';




//////////////////////////////
// lens options
$html .= '<div id="lens_div" class="tab-content current">';
$json = json_decode(file_get_contents('form_json/options_lens.json'),true);

$html .= '<div class="slot">';
$html .= '<p class="title">Lenses</p>';
$html .= createInputs($json['physical']['input_list'],$json['physical']['type'],'physical');

$html .= '<hr class="myhr">';
$html .= '<div id="lens-1" class="lens">';
$html .= '<h4>Lens 1</h4>';
$html .= '<button type="button" class="remove_lens">Remove</button>';
$html .= createSelect($json['lenses']['select_list'],'lenses');
$html .= '<hr class="myhr">';
$html .= '</div>';
$html .= '<button type="button" id="add_lens">Add lens</button>';
$html .= '</div>';

$html .= '<div class="slot">';
$html .= '<p class="title">Perturbations</p>';
$html .= createSelect($json['perturbations']['select_list'],'perturbations');
$html .= '</div>';

$html .= '<div class="slot">';
$html .= '<p class="title">Generic options</p>';
$html .= createInputs($json['print_all']['input_list'],$json['print_all']['type'],'print_all');
$html .= '</div>';

$html .= '</div>';




//////////////////////////////
// observation options
$html .= '<div id="observer_div" class="tab-content">';
$json = json_decode(file_get_contents('form_json/options_observer.json'),true);

$html .= '<div class="slot">';
$html .= '<p class="title">Field of view</p>';
$html .= createInputs($json['iplane']['input_list'],$json['iplane']['type'],'iplane');
$html .= '</div>';

$html .= '<div class="slot">';
$html .= '<p class="title">PSF</p>';
$html .= createInputs($json['psf']['input_list'],$json['psf']['type'],'psf');
$html .= '</div>';

$html .= '<div class="slot">';
$html .= '<p class="title">Noise</p>';
$html .= createSelect($json['noise']['select_list'],'noise');
$html .= '</div>';

$html .= '<div class="slot">';
$html .= '<p class="title">Mask</p>';
$html .= createSelect($json['mask']['select_list'],'mask');
$html .= '</div>';

$html .= '</div>';

$html .= '</div>';
$html .= '</form>';


$html .= '<button type="button" id="create">Create!</button>';
$html .= '<pre id="json"></pre>';


$html .='</body>';
$html .='</html>';

file_put_contents('index.html',$html);
//echo $html;






function createSelect($arr,$id){
    $str = '<select class="myselect" id="'.$id.'">';
    $str .= '<option disabled selected="selected">Select option</option>';
    for($i=0;$i<count($arr);$i++){
	$str .= '<option value="'.$arr[$i]['name'].'">'.$arr[$i]['label'].'</option>';
    }
    $str .= '</select>';

    for($i=0;$i<count($arr);$i++){
	$str .= '<div class="myselect_option '.$arr[$i]['name'].'" >';
	$str .= '<p class="legend">' . $arr[$i]['description'] . '</p>';
	$str .= createInputs($arr[$i]['input_list'],$arr[$i]['type'],$arr[$i]['name']);
	$str .= '</div>';
    }

    return $str;
}


function createInputs($arr,$type,$id){
    if( $type == 'nlpar_list' ){

	$str = '<table class="input_table nlpars" id="'.$id.'">';
	//$str .= '<thead><th>&nbsp;</th><th>Min.</th><th>Value</th><th>Max.</th><th>Fix</th><th>&nbsp;</th></thead>';
	$str .= '<thead><th>&nbsp;</th><th>&nbsp;</th><th>Value</th><th>&nbsp;</th></thead>';
	for($i=0;$i<count($arr);$i++){
	    $str .= '<tr>';
	    $str .= '<td><label>'.$arr[$i]['label'].'</label></td>';
	    $str .= '<td><input type="hidden" name="'.$arr[$i]['name'].'_nam" class="nlpar_value" value="'.$arr[$i]['name'].'" /></td>';
	    //$str .= '<td><input type="number" step="any" name="'.$arr[$i]['name'].'_min" class="nlpar_range" placeholder="Min" /></td>';
	    $str .= '<td><input type="number" step="any" name="'.$arr[$i]['name'].'_val" class="nlpar_value" placeholder="Value" /></td>';
	    //$str .= '<td><input type="number" step="any" name="'.$arr[$i]['name'].'_max" class="nlpar_range" placeholder="Max" /></td>';
	    //$str .= '<td><input type="checkbox" name="'.$arr[$i]['name'].'_fix" class="nlpar_fix" /></td>';
	    $str .= '<td><p class="legend">'.$arr[$i]['legend'].'</p></td>';
	    $str .= '</tr>';
	}
	$str .= '</table>';

    } else {

	$str = '<table class="input_table" id="'.$id.'">';
	for($i=0;$i<count($arr);$i++){
	    $label  = $arr[$i]['label'];
	    $legend = $arr[$i]['legend'];
	    unset($arr[$i]['label']);
	    unset($arr[$i]['legend']);

	    if( $arr[$i]['type'] != "hidden" ){
		$str .= '<tr>';
		$str .= '<td><label>'.$label.'</label></td>';
		if( $arr[$i]['type'] == 'select' ){
		    
		} else {
		    $str .= '<td><input ';
		    foreach($arr[$i] as $key => $val){
			$str .= $key . '="' .$val. '" ';
		    }
		    $str .= ' /></td>';
		}
		$str .= '<td><p class="legend">'.$legend.'</p></td>';
		$str .= '</tr>';
	    } else {
		$str .= '<tr><td></td><td><input ';
		foreach($arr[$i] as $key => $val){
		    $str .= $key . '="' .$val. '" ';
		}
		$str .= ' /></td><td></td></tr>';
	    }
	}
	$str .= '</table>';

    }
    return $str;
}




/*
   function createContent($json_file){
   $json = json_decode(file_get_contents($json_file),true);
   $cnt = '';
   foreach($json as $key => $obj){
   if( $json[$key]["type"] == "select" ){
   $cnt .= '<div class="myselect">';
   foreach($json[$key]["select_list"] as $key2 => $obj2){
   $cnt .= $obj2["name"] . "," ;
   }
   $cnt .= "</div>";
   } else {
   $cnt .= '<div class="myinput">';
   $cnt .= $key . '<br>';
   foreach($json[$key]["input_list"] as $key2 => $obj2){
   $cnt .= $obj2["name"] . "," ;
   }	    
   $cnt .= "</div>";
   }
   }


   return $cnt;
   }
 */
?>
