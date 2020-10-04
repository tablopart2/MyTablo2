use hokimdb

create table lotto_database
(
	num decimal(9,0) 
	,publishdate datetime not null
	,numwin1 decimal(16,0) 
	,moneywin1 decimal(16,0) 
	,numwin2 decimal(16,0) 
	,moneywin2 decimal(16,0) 
	,numwin3 decimal(16,0) 
	,moneywin3 decimal(16,0) 
	,numwin4 decimal(16,0) 
	,moneywin4 decimal(16,0) 
	,numwin5 decimal(16,0) 
	,moneywin5 decimal(16,0) 

	,first_n decimal(16,0) 
	,second_n decimal(16,0) 
	,third_n decimal(16,0) 
	,fourth_n decimal(16,0) 
	,fifth_n decimal(16,0) 
	,sixth_n decimal(16,0) 
	,bonus_n decimal(16,0) 
	primary key(num, publishdate)
)

