% % clear previous settings
% setdbprefs('JDBCDataSourceFile', '')

% Database Server
host = 'TUD276632:3306'; % 'localhost:3333';

% Database Username/Password
user = 'abun';
password = 'welcome';

% Database Name
dbName = 'landfills2'; 

% JDBC Parameters
jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
jdbcDriver = 'com.mysql.jdbc.Driver';

% Set this to the path to your MySQL Connector/J JAR
javaaddpath('./mysql-connector-java-5.1.25/mysql-connector-java-5.1.25-bin.jar')

% Create the database connection object
dbConn = database(dbName, user, password, jdbcDriver, jdbcString);
isconnection(dbConn);

% add paths of subsequent folders to workspace
addpath(genpath('../'));