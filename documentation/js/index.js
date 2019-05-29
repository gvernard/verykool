$(document).ready(function(){
    //$("#tabs").tabs();

    $('ul.tabs li').click(function(){
	var tab_id = $(this).attr('data-tab');
	
	$('ul.tabs li').removeClass('current');
	$('.tab-content').removeClass('current');
	
	$(this).addClass('current');
	$("#"+tab_id).addClass('current');
    })


    $('.myselect_option').hide();
    $('.myselect').find('option:eq(0)').prop('selected',true);

    $('body').on('change','.myselect',function(){
	$(this).siblings('.selected').hide().removeClass('selected');
	$(this).siblings('.'+$(this).val()).show().addClass('selected');
    });


    var $lens = $('#lens-1').clone();
    $('#add_lens').click(function(){
	var $new_lens = $lens.clone();
	$(this).before($new_lens);
	renumber_lenses();
    });

    $('#lens_div').on('click','.remove_lens',function(){
	$(this).parent().remove();
	renumber_lenses();
    });

    var $anal_source = $('#anal-source-1').clone();
    $('#add_source').click(function(){
	var $new_source = $anal_source.clone();
	$(this).before($new_source);
	renumber_sources();
    });

    $('#source_div').on('click','.remove_source',function(){
	$(this).parent().remove();
	renumber_sources();
    });



    $('#create').click(function(){
	var myjson = createJSON();
	$('#json').text( JSON.stringify(myjson,undefined,2) );
    });
});

function renumber_lenses(){
    $('.lens').each(function(index,element){
	$(element).attr('id','lens-'+index);
	$(element).find('h4').text('Lens '+(index+1));
    });
}

function renumber_sources(){
    $('.anal-source').each(function(index,element){
	$(element).attr('id','source-'+index);
	$(element).find('h4').text('Analytic Source '+(index+1));
    });
}

function get_inputs(id){
    var obj = {};
    $('#'+id).find('input').each(function(){
	var name = $(this).attr('name');
	var val  = $(this).val();
        if( !isNaN(val) ){
            val = parseFloat(val);
        }
	obj[name] = val;
    });
    return obj;
};

function get_nlpars(id){
    var arr = [];
    $('#'+id).find('tbody > tr').each(function(){
	var obj = {};
	var name = $(this).find('input:hidden').val();
	obj['nam'] = name;
	var val  = parseFloat( $(this).find('input[name="'+name+'_val"]').val() );
	obj['val'] = val;
	arr.push( obj );
    });
    return arr;
}


function createJSON(){
    var myjson = {};
    
    // Image plane
    var iplane = get_inputs('iplane');
    myjson['general'] = {};
    myjson['general']['iplane'] = iplane;
    
    // Perturbations
    
    // FProject
    var fproject = {};
    //var print_all = get_inputs('print_all');
    fproject['physical'] = get_nlpars('physical');
    
    var sel_source = $('#source').val();
    var source = get_inputs(sel_source);
    source['type'] = sel_source;
    fproject['source'] = source;
    
    
    myjson['fproject'] = fproject;
    
    
    
    // Noise
    var sel_noise = $('#noise').val();
    var noise = get_inputs(sel_noise);
    noise['noise_flag'] = sel_noise;
    myjson['noise'] = noise;
    
    // PSF
    var psf = get_inputs('psf');
    myjson['addmachine'] = psf;
    
    // Mask
    var sel_mask = $('#mask').val();
    var mask = get_inputs(sel_mask);
    mask['mask_flag'] = sel_mask;
    myjson['automask'] = mask;

    return myjson;
}
