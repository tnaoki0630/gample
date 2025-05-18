#import "Init.h"
#import "Constant.h"
#import <string>
#import <fstream>
#import <iostream>
#import <sstream>
#import <vector>

@interface Init() {
    // Private properties to store parsed values
    struct FragForEquation _fragEquation;
    struct ParamForTimeIntegration _timeIntegration;
    std::vector<struct ParamForParticle> _particles;
    std::vector<struct BoundaryConditionForParticle> _particleBoundaries;
    std::vector<struct SourceForParticle> _particleSources;
    struct ParamForField _field;
    std::vector<struct BoundaryConditionForField> _fieldBoundaries;
    // optional values
    std::vector<bool> _flagForPNum;
    std::vector<int> _initPNumPerCell;
    std::vector<float> _initWeightFromDens;
    std::vector<float> _currentDensity;
}

@property (nonatomic, strong) NSString *filePath;

@end

@implementation Init

- (instancetype)parseInputFile:(NSString*)InputFilePath {
    self = [super init];
    if (self) {
        _filePath = InputFilePath;
        
        // Initialize vectors
        _particles = std::vector<struct ParamForParticle>();
        _particleBoundaries = std::vector<struct BoundaryConditionForParticle>();
        _particleSources = std::vector<struct SourceForParticle>();
        _fieldBoundaries = std::vector<struct BoundaryConditionForField>();
        _flagForPNum = std::vector<bool>();
        _initPNumPerCell = std::vector<int>();
        _initWeightFromDens = std::vector<float>();
        _currentDensity = std::vector<float>();

        // default param for field
        _field.ngx = -1;
        _field.ngy = -1;
        _field.dx = -1.0;
        _field.dy = -1.0;
        _field.InitTypeE = NULL;
        _field.InitTypeB = NULL;
        _field.ampE[0] = 0.0;
        _field.ampE[1] = 0.0;
        _field.ampE[2] = 0.0;
        _field.ampB[0] = 0.0;
        _field.ampB[1] = 0.0;
        _field.ampB[2] = 0.0;
        _field.FilePathEx = NULL;
        _field.FilePathEy = NULL;
        _field.FilePathEz = NULL;
        _field.FilePathBx = NULL;
        _field.FilePathBy = NULL;
        _field.FilePathBz = NULL;
        _field.weightOrder = -1;
        _field.ngb = -1;
        _field.maxiter = 200;
        _field.tolerance = 1e-5;
        
        // Convert NSString to std::string for C++ file handling
        std::string filePath = [InputFilePath UTF8String];
        std::ifstream inputFile(filePath);
        
        if (!inputFile.is_open()) {
            NSLog(@"Error: Could not open input file %@", InputFilePath);
            return nil;
        }
        
        // Parse the file
        if (![self parseFile:inputFile]) {
            NSLog(@"Error: Failed to parse input file %@", InputFilePath);
            return nil;
        }
        
        inputFile.close();
    }
    return self;
}

- (BOOL)parseFile:(std::ifstream&)inputFile {
    std::string line;
    std::string currentSection = "";
    
    while (std::getline(inputFile, line)) {
        // Skip empty lines
        if (line.empty()) continue;
        
        // Skip comment lines
        if (line[0] == '#') continue;
        
        // Check if this is a section header
        if (line == "FlagForEquation") {
            currentSection = line;
            continue;
        } else if (line == "ParamForTimeIntegration") {
            currentSection = line;
            continue;
        } else if (line == "ParamForParticle") {
            currentSection = line;
            continue;
        } else if (line == "BoundaryConditionForParticle") {
            currentSection = line;
            continue;
        } else if (line == "SourceForParticle") {
            currentSection = line;
            continue;
        } else if (line == "ParamForField") {
            currentSection = line;
            continue;
        } else if (line == "BoundaryConditionForField") {
            currentSection = line;
            continue;
        } else if (line == "GOGO") {
            // End of file
            break;
        }
        
        // Parse based on current section
        if (currentSection == "FlagForEquation") {
            [self parseFlagForEquation:line];
        } else if (currentSection == "ParamForTimeIntegration") {
            [self parseParamForTimeIntegration:line];
        } else if (currentSection == "ParamForParticle") {
            if (line == "/") {
                // End of a particle definition
                continue;
            }
            
            [self parseParamForParticle:line inputFile:inputFile];
        } else if (currentSection == "BoundaryConditionForParticle") {
            if (line == "/") {
                // End of a boundary condition
                continue;
            }
            [self parseBoundaryConditionForParticle:line inputFile:inputFile];
        } else if (currentSection == "SourceForParticle") {
            if (line == "/") {
                // End of a source condition
                continue;
            }
            [self parseSourceForParticle:line inputFile:inputFile];
        } else if (currentSection == "ParamForField") {
            [self parseParamForField:line];
        } else if (currentSection == "BoundaryConditionForField") {
            if (line == "/") {
                // End of a boundary condition
                continue;
            }
            [self parseBoundaryConditionForField:line inputFile:inputFile];
        }
    }
    
    return YES;
}

- (void)parseFlagForEquation:(const std::string&)line {
    std::istringstream iss(line);
    std::string key;
    int value;
    
    if (iss >> key >> value) {
        if (key == "Particle") {
            // number of specimens
            _fragEquation.Particle = value;
        } else if (key == "EMField") {
            // (developing)0: invariant EMField
            // 1: Electro-static field w/ poisson eq.
            // (developing)2: Electro-Magnetic field w/ FDTD
            _fragEquation.EMField = value;
        } else if (key == "MCCollision") {
            // 0: collisionless
            // (developing)1: MCC-collision
            // (developing)2: null-collision
            _fragEquation.MCCollision = value;
        }
    }
}

- (void)parseParamForTimeIntegration:(const std::string&)line {
    std::istringstream iss(line);
    std::string key;
    double value;
    
    if (iss >> key >> value) {
        if (key == "End") {
            _timeIntegration.EndCycle = (int)value;
        } else if (key == "ParticleOutput") {
            _timeIntegration.ptclOutCycle = (int)value;
        } else if (key == "FieldOutput") {
            _timeIntegration.fldOutCycle = (int)value;
        } else if (key == "TimeStep") {
            _timeIntegration.dt = value;
        }
    }
}

- (void)parseParamForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile {
    // Start a new particle definition
    struct ParamForParticle particle;
    // Initialize fields to default values
    particle.pName = NULL;
    particle.pNum = 0;
    particle.pNumMax = 0;
    particle.q = 0.0;
    particle.m = 0.0;
    particle.w = 0.0;
    particle.genType = NULL;
    particle.genX[0] = 0.0;
    particle.genX[1] = 0.0;
    particle.genY[0] = 0.0;
    particle.genY[1] = 0.0;
    particle.genU[0] = 0.0;
    particle.genU[1] = 0.0;
    particle.genU[2] = 0.0;
    particle.genT = 0.0;
    _flagForPNum.push_back(false);
    _initPNumPerCell.push_back(0);
    
    std::istringstream iss(line);
    std::string key;
    std::string value;
    if (iss >> key >> value) {
        if (key == "ParticleName") {
            particle.pName = [NSString stringWithUTF8String:value.c_str()];
        }
    }
    

    // Parse subsequent lines for this particle until we hit a '/' line
    std::string valueLine;
    while (std::getline(inputFile, valueLine)) {
        if (valueLine == "/") {
            // End of particle definition
            break;
        }
        
        std::istringstream lineIss(valueLine);
        std::string paramKey;
        std::string paramValue;
        std::string val1;
        std::string val2;
        
        if (lineIss >> paramKey >> paramValue) {
            // NSLog(@"paramKey = %s, paramValue = %s", paramKey.c_str(), paramValue.c_str());
            if (paramKey == "InitialParticleNumber") {
                particle.pNum = std::stoi(paramValue);
                _flagForPNum.back() = true;
            } else if (paramKey == "InitPtclNumPerCell") {
                // 初期の particlePerCell で指定
                _initPNumPerCell.back() = std::stoi(paramValue);
            } else if (paramKey == "MaxParticleNumber") {
                particle.pNumMax = std::stoi(paramValue);
            } else if (paramKey == "Charge") {
                particle.q = std::stod(paramValue);
            } else if (paramKey == "Mass") {
                particle.m = std::stod(paramValue);
            } else if (paramKey == "Weight[1/cm]") {
                particle.w = std::stod(paramValue);
                _initWeightFromDens.push_back(0);
            } else if (paramKey == "WeightFromDens[1/cm3]") {
                // 初期の数密度で指定
                _initWeightFromDens.push_back(std::stod(paramValue));
            } else if (paramKey == "GenerateType") {
                particle.genType = [NSString stringWithUTF8String:paramValue.c_str()];
            } else if (paramKey == "InitialPosX") {
                if (paramValue == "auto"){
                    particle.genX[0] = -1.0;
                    particle.genX[1] = -1.0;
                }else{
                    particle.genX[0] = std::stod(paramValue);
                    if (lineIss >> val1) {
                        particle.genX[1] = std::stod(val1);
                    }
                }
            } else if (paramKey == "InitialPosY") {
                if (paramValue == "auto"){
                    particle.genY[0] = -1.0;
                    particle.genY[1] = -1.0;
                }else{
                    particle.genY[0] = std::stod(paramValue);
                    if (lineIss >> val1) {
                        particle.genY[1] = std::stod(val1);
                    }
                }
            } else if (paramKey == "InitialVel") {
                particle.genU[0] = std::stod(paramValue);
                if (lineIss >> val1 >> val2) {
                    particle.genU[1] = std::stod(val1);
                    particle.genU[2] = std::stod(val2);
                }
            } else if (paramKey == "InitialTemp[eV]") {
                particle.genT = std::stod(paramValue)*evtok;
            }
        }
    }    
    // Add the completed particle to vector
    _particles.push_back(particle);
}

- (void)parseBoundaryConditionForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile {
    struct BoundaryConditionForParticle boundary;
    // Initialize fields to default values
    boundary.position = NULL;
    boundary.type = NULL;
    boundary.val = 0.0;
    
    // Parse the first line (RegionName)
    std::istringstream nameIss(line);
    std::string key, value;
    if (nameIss >> key >> value && key == "RegionName") {
        boundary.position = [NSString stringWithUTF8String:value.c_str()];;
    }
    
    // Parse the type line
    std::string valueLine;
    if (std::getline(inputFile, valueLine)) {
        std::istringstream typeIss(valueLine);
        if (typeIss >> key >> value && key == "BCType") {
            boundary.type = [NSString stringWithUTF8String:value.c_str()];;
        }
    }
    
    // Parse optional value line (if any before the '/' delimiter)
    if (std::getline(inputFile, valueLine) && valueLine != "/") {
        std::istringstream valueIss(valueLine);
        double val;
        if (valueIss >> val) {
            boundary.val = val;
        }
        // Skip the next line if it's a '/'
        std::string nextLine;
        if (std::getline(inputFile, nextLine) && nextLine != "/") {
            // Put back if not a delimiter
            inputFile.seekg(-nextLine.length() - 1, std::ios_base::cur);
        }
    }
    
    // Add the boundary to vector
    _particleBoundaries.push_back(boundary);
}

- (void)parseSourceForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile {
    struct SourceForParticle source;
    // Initialize fields to default values
    source.pName = NULL;
    source.genType = NULL;
    source.src = 0.0;
    source.genX[0] = 0.0;
    source.genX[1] = 0.0;
    source.genY[0] = 0.0;
    source.genY[1] = 0.0;
    source.genU[0] = 0.0;
    source.genU[1] = 0.0;
    source.genU[2] = 0.0;
    source.genT = 0.0;
   
    std::istringstream iss(line);
    std::string key;
    std::string value;
    if (iss >> key >> value) {
        if (key == "ParticleName") {
            source.pName = [NSString stringWithUTF8String:value.c_str()];
        }
    }

    std::string valueLine;
    while (std::getline(inputFile, valueLine)) {
        if (valueLine == "/") {
            // End of particle definition
            break;
        }
        
        std::istringstream lineIss(valueLine);
        std::string paramKey;
        std::string paramValue;
        std::string val1;
        std::string val2;
        
        // NSLog(@"ParamValue: %@", [NSString stringWithUTF8String:line.c_str()]);
        if (lineIss >> paramKey >> paramValue) {
            if (paramKey == "GenerateType") {
                source.genType = [NSString stringWithUTF8String:paramValue.c_str()];
                if ([source.genType isEqualToString:@"hollow-cathode"]){
                    _currentDensity.push_back(0); // 使わないので0埋め
                }
            } else if (paramKey == "SourceValue[1/(cm s)]") {
                source.src = std::stod(paramValue);
                _currentDensity.push_back(0); // 使わないので0埋め
            } else if (paramKey == "SourceValue[A/m2]") {
                _currentDensity.push_back(std::stod(paramValue)*JtosJ);
            } else if (paramKey == "GeneratePosX") {
                if (paramValue == "auto"){
                    source.genX[0] = -1.0;
                    source.genX[1] = -1.0;
                }else{
                    source.genX[0] = std::stod(paramValue);
                    if (lineIss >> val1) {
                        source.genX[1] = std::stod(val1);
                    }
                }
            } else if (paramKey == "GeneratePosY") {
                if (paramValue == "auto"){
                    source.genY[0] = -1.0;
                    source.genY[1] = -1.0;
                }else{
                    source.genY[0] = std::stod(paramValue);
                    if (lineIss >> val1) {
                        source.genY[1] = std::stod(val1);
                    }
                }
            } else if (paramKey == "GenerateVel") {
                source.genU[0] = std::stod(paramValue);
                if (lineIss >> val1 >> val2) {
                    source.genU[1] = std::stod(val1);
                    source.genU[2] = std::stod(val2);
                }
            } else if (paramKey == "GenerateTemp[eV]") {
                source.genT = std::stod(paramValue)*evtok;
            }
        }
    }    
    // Add the completed source to vector
    _particleSources.push_back(source);
}

- (void)parseParamForField:(const std::string&)line {
    
    std::istringstream iss(line);
    std::string key;
    std::string value;
    std::string val1;
    std::string val2;
    
    if (iss >> key >> value) {
        if (key == "NumberOfGridX") {
            _field.ngx = stoi(value);
        } else if (key == "NumberOfGridY") {
            _field.ngy = stoi(value);
        } else if (key == "GridSizeOfX") {
            _field.dx = stod(value);
        } else if (key == "GridSizeOfY") {
            _field.dy = stod(value);
        } else if (key == "InitializeTypeOfE") {
            _field.InitTypeE = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "AmplitudeOfE[V/m]") {
            _field.ampE[0] = stod(value)*VtoG;
            if (iss >> val1 >> val2) {
                _field.ampE[1] = stod(val1)*VtoG;
                _field.ampE[2] = stod(val2)*VtoG;
            }
        } else if (key == "FilePathOfEx") {
            _field.FilePathEx = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "FilePathOfEy") {
            _field.FilePathEy = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "FilePathOfEz") {
            _field.FilePathEz = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "InitializeTypeOfB") {
            _field.InitTypeB = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "AmplitudeOfB[V/m]") {
            _field.ampB[0] = stod(value)*TtoG;
            if (iss >> val1 >> val2) {
                _field.ampB[1] = stod(val1)*TtoG;
                _field.ampB[2] = stod(val2)*TtoG;
            }
        } else if (key == "FilePathOfBx") {
            _field.FilePathBx = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "FilePathOfBy") {
            _field.FilePathBy = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "FilePathOfBz") {
            _field.FilePathBz = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "WeightingOrder") {
            _field.weightOrder = stoi(value);
            _field.ngb = stoi(value)/2; // 5th-order -> ngb = 2
        } else if (key == "MaxIterForPoisson") {
            _field.maxiter = stoi(value);
        } else if (key == "TolForPoisson") {
            _field.tolerance = stof(value);
        }
    }
}

- (void)parseBoundaryConditionForField:(const std::string&)line inputFile:(std::ifstream&)inputFile {
    struct BoundaryConditionForField boundary;
    // Initialize fields to default values
    boundary.position = NULL;
    boundary.type = NULL;
    boundary.val = 0.0;
    
    // Parse the first line (RegionName)
    std::istringstream nameIss(line);
    std::string key, value;
    if (nameIss >> key >> value && key == "RegionName") {
        boundary.position = [NSString stringWithUTF8String:value.c_str()];
    }
    
    // Parse the type line
    std::string valueLine;
    if (std::getline(inputFile, valueLine)) {
        std::istringstream typeIss(valueLine);
        if (typeIss >> key >> value && key == "BCType") {
            boundary.type = [NSString stringWithUTF8String:value.c_str()];
        }
    }
    
    // Parse optional value line (if any before the '/' delimiter)
    if (std::getline(inputFile, valueLine) && valueLine != "/") {
        std::istringstream valueIss(valueLine);
        double val;
        if (valueIss >> val) {
            boundary.val = val*VtosV;
        }
        // Skip the next line if it's a '/'
        std::string nextLine;
        if (std::getline(inputFile, nextLine) && nextLine != "/") {
            // Put back if not a delimiter
            inputFile.seekg(-nextLine.length() - 1, std::ios_base::cur);
        }
    }
    
    // Add the boundary to vector
    _fieldBoundaries.push_back(boundary);
}


- (BOOL)checkInput{
    // check fieldParam
    if (_field.ngx > 0 && _field.ngy > 0 && _field.ngb >= 0 && _field.dx > 0.0 && _field.dy > 0.0){
        // check vectorsize
        if(_particles.size() != _initPNumPerCell.size()){
            NSLog(@"check vectorsize failed.");
            NSLog(@"_particles.size: %zu", _particles.size());
            NSLog(@"_initPNumPerCell.size: %zu", _initPNumPerCell.size());
            return false;
        }
        if(_particleSources.size() != _currentDensity.size()){
            NSLog(@"check vectorsize failed.");
            NSLog(@"_particleSources.size: %zu", _particleSources.size());
            NSLog(@"_currentDensity.size: %zu", _currentDensity.size());
            return false;
        }

        // check particleParam
        // iteration
        for (int s = 0; s < _particles.size(); s++){
            // pNumMax は初期化必須
            if(_particles[s].pNumMax == 0){
                NSLog(@"check particleParam for %@ failed.", _particles[s].pName);
                return false;
            }
            // pNum か initPNumPerCell のどちらかは初期化必須
            if(_flagForPNum[s]){
                // OK
            }else if(_initPNumPerCell[s] > 0){
                // culculate pNum from initPNumPerCell
                _particles[s].pNum = _initPNumPerCell[s]*_field.ngx*_field.ngy;
            }else{
                NSLog(@"check particleParam for %@ failed.", _particles[s].pName);
                return false;
            }
            NSLog(@"check pNum: %d", _particles[s].pNum);
            // weight か _initWeightFromDens のどちらかは初期化必須
            if(_particles[s].w > 0.0){
                // OK
            }else if(_initWeightFromDens[s] > 0.0){
                // culculate weight from initial number density
                float ppc = (float)_particles[s].pNum/(float)(_field.ngx*_field.ngy);
                _particles[s].w = _initWeightFromDens[s]*_field.dx*_field.dy/ppc;
            }else{
                NSLog(@"check particleParam for %@ failed.", _particles[s].pName);
                return false;
            }

            // check Source
            // iteration
            for (size_t i = 0; i < _particleSources.size(); i++){
                // 粒子種が一致し、かつ genType が hollow-cathode でない生成条件をチェック
                if([_particles[s].pName isEqualToString:_particleSources[i].pName] && ![_particleSources[i].genType isEqualToString:@"hollow-cathode"]){
                    if(_particleSources[i].src > 0.0){
                        // OK
                    }else if(_currentDensity[i] > 0.0){
                        // culculate src from currentDensity
                        _particleSources[i].src = _currentDensity[i]*_field.ngy*_field.dy/(ec*_particles[s].w);
                    }else{
                        NSLog(@"check particleSources for %@ failed.", _particles[s].pName);
                        return false;
                    }
                }
            }
        }
    }else{
        NSLog(@"check fieldparam failed.");
        return false;
    }
    return true;
}


// Public accessor methods to get the parsed data
- (struct FragForEquation)getFragForEquation {
    return _fragEquation;
}

- (struct ParamForTimeIntegration)getParamForTimeIntegration {
    return _timeIntegration;
}

- (NSArray*)getParamForParticle {
    NSMutableArray *result = [NSMutableArray array];
    for (const auto& particle : _particles) {
        // Create a copy to return
        struct ParamForParticle copy = particle;
        [result addObject:[NSValue value:&copy withObjCType:@encode(struct ParamForParticle)]];
    }
    return result;
}

- (NSArray*)getParticleBoundaries {
    NSMutableArray *result = [NSMutableArray array];
    for (const auto& boundary : _particleBoundaries) {
        // Create a copy to return
        struct BoundaryConditionForParticle copy = boundary;
        [result addObject:[NSValue value:&copy withObjCType:@encode(struct BoundaryConditionForParticle)]];
    }
    return result;
}

- (NSArray*)getParticleSources {
    NSMutableArray *result = [NSMutableArray array];
    for (const auto& source : _particleSources) {
        // Create a copy to return
        struct SourceForParticle copy = source;
        [result addObject:[NSValue value:&copy withObjCType:@encode(struct SourceForParticle)]];
    }
    return result;
}

- (struct ParamForField)getParamForField {
    return _field;
}

- (NSArray*)getFieldBoundaries {
    NSMutableArray *result = [NSMutableArray array];
    for (const auto& boundary : _fieldBoundaries) {
        // Create a copy to return
        struct BoundaryConditionForField copy = boundary;
        [result addObject:[NSValue value:&copy withObjCType:@encode(struct BoundaryConditionForField)]];
    }
    return result;
}
@end