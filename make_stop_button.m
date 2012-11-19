function make_stop_button()

  % Initialize variables and graphics:

  keepLooping = true;
  hFigure = figure;
  hButton = uicontrol('Style','pushbutton','Parent',hFigure,...
                      'String','Stop','Callback',@stop_fcn);

  % Keep looping until button is pressed:

  while keepLooping,
    drawnow;
  end

  % Delete figure:

  delete(hFigure);

  % Nested callback:

  function stop_fcn(source,event)
    keepLooping = false;
  end

end