function Group_Migrate(Community1,Community2,migrate_num)
    exchange_set = randi(min(Community1.num,Community2.num),migrate_num);
    Community1.individual{exchange_set} = copy(Community2.individual...
        {exchange_set});
end