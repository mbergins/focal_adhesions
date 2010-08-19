function cyto_set_check(handles, nr)

if nr == 1
    set(handles.correlate(1),'Check','on');
else
    set(handles.correlate(1),'Check','off'); 
end
if nr == 2
    set(handles.subtract_v(1),'Check','on');
else
    set(handles.subtract_v(1),'Check','off');
end
if nr == 3
    set(handles.inverse(1),'Check','on');
else
    set(handles.inverse(1),'Check','off');
end
if nr == 4
    set(handles.subtract(1),'Check','on');
else
    set(handles.subtract(1),'Check','off');
end
if nr == 5
    set(handles.multiply(1),'Check','on');
else
    set(handles.multiply(1),'Check','off');    
end
if nr == 6
    set(handles.divide(1),'Check','on');
else
    set(handles.divide(1),'Check','off');
end
if nr == 7
    set(handles.get_bleach(1),'Check','on');
else
    set(handles.get_bleach(1),'Check','off');
end
if nr == 8
    set(handles.img_vec_mult(1),'Check','on');
else
    set(handles.img_vec_mult(1),'Check','off');
end
if nr == 9
    set(handles.rac(1),'Check','on');
else
    set(handles.rac(1),'Check','off');
end
if nr == 10
    set(handles.activity_from_edge(1),'Check','on');
else
    set(handles.activity_from_edge(1),'Check','off');
end
if nr == 11
    set(handles.automatic(1),'Check','on');
else
    set(handles.automatic(1),'Check','off');
end
if nr == 12
    set(handles.manual(1),'Check','on');
else
    set(handles.manual(1),'Check','off');
end
if nr == 13
    set(handles.bleach_correct(1),'Check','on');
else
    set(handles.bleach_correct(1),'Check','off');
end
if nr == 14
    set(handles.curvature(1),'Check','on');
else
    set(handles.curvature(1),'Check','off');
end
if nr == 15
    set(handles.subtract_bg(1),'Check','on');
else
    set(handles.subtract_bg(1),'Check','off');
end
if nr == 16
    set(handles.translate(1),'Check','on');
else
    set(handles.translate(1),'Check','off');
end
if nr == 17
    set(handles.get_translation(1),'Check','on');
else
    set(handles.get_translation(1),'Check','off');
end
if nr == 18
    set(handles.make_movie(1),'Check','on');
else
    set(handles.make_movie(1),'Check','off');
end
if nr == 19
    set(handles.average_image(1),'Check','on');
else
    set(handles.average_image(1),'Check','off');
end
if nr == 20
    set(handles.overlay(1),'Check','on');
else
    set(handles.overlay(1),'Check','off');
end
if nr == 21
    set(handles.average_vector_field(1),'Check','on');
else
    set(handles.average_vector_field(1),'Check','off');
end
if nr == 22
    set(handles.mpm_to_vec(1),'Check','on');
else
    set(handles.mpm_to_vec(1),'Check','off');
end
if nr == 23
    set(handles.track_analysis(1),'Check','on');
else
    set(handles.track_analysis(1),'Check','off');
end
if nr == 24
    set(handles.scores_to_vec(1),'Check','on');
else
    set(handles.scores_to_vec(1),'Check','off');
end
if nr == 25
    set(handles.analyse_prot(1),'Check','on');
else
    set(handles.analyse_prot(1),'Check','off');
end
if nr == 26
    set(handles.pka(1),'Check','on');
else
    set(handles.pka(1),'Check','off');
end