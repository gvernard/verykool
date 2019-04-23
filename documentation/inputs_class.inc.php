<?php 
class myInput{
    var $name;
    var $type;
    var $value;
    var $def_value;
    var $descr;
    
    function myInput($name,$address,$age){
	$this->name = $name;
	$this->address = $address;
	$this->age = $age;
    }

    function printInput(){
	$str = '<tr>';
	$str .= '<td><label></label></td>';
	$str .= '<td><input name="" type="" value="" /></td>';
	$str .= '<td><span class="description">'.$this->descr.'</span></td>';
	$str .= '</tr>';
    }
}


?>
