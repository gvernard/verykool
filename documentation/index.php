<?php 

$html ='<html>';
$html .='<head>';
$html .='</head>';

$html .='<body>';

$html .='<h2>Tool to generate input file for mocks</h2>';


// General options
$html .= '<div id="general">';
$html .= '<h4>General</h4>';
$html .= '</div>';


// Perturbations options
$html .= '<div id="perturbations">';
$html .= '<h4>Perturbations</h4>';
$html .= '</div>';


// FProject options
$html .= '<div id="fproject">';
$html .= '<h4>FProject</h4>';
$html .= '</div>';



$html .='</body>';
$html .='</html>';

file_put_contents('index.html',$html);
//echo $html;
?>
