function logging_func(message_str)

    %
    % logging_func
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %

    str_datetime = string(datetime('now','Format','MM/dd HH:mm:ss.SSS'));
    fprintf("[%s] %s\n",str_datetime,message_str);
end