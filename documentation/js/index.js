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
    $('.myselect').change(function(){
	$(this).siblings('.selected').hide().removeClass('selected');
	$(this).siblings('.'+$(this).val()).show().addClass('selected');
    });

});
