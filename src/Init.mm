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
    struct ParamForField _field;
    std::vector<struct BoundaryConditionForField> _fieldBoundaries;
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
        _fieldBoundaries = std::vector<struct BoundaryConditionForField>();
        
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
            _fragEquation.Particle = value;
        } else if (key == "EMField") {
            _fragEquation.EMField = value;
        } else if (key == "MCCollision") {
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
            _timeIntegration.pOutCycle = (int)value;
        } else if (key == "FieldOutput") {
            _timeIntegration.fOutCycle = (int)value;
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
    particle.GenerateType = NULL;
    particle.initX[0] = 0.0;
    particle.initX[1] = 0.0;
    particle.initY[0] = 0.0;
    particle.initY[1] = 0.0;
    particle.initU[0] = 0.0;
    particle.initU[1] = 0.0;
    particle.initU[2] = 0.0;
    particle.initT = 0.0;
    
    std::istringstream iss(line);
    std::string key;
    std::string value;
    if (iss >> key >> value) {
        if (key == "ParticleName") {
            particle.pName = [NSString stringWithUTF8String:value.c_str()];
            NSLog(@"value: %@", [NSString stringWithUTF8String:value.c_str()]);
            NSLog(@"pName: %@", particle.pName);
        }
    }
    
    // Parse subsequent lines for this particle until we hit a '/' line
    std::string particleLine;
    while (std::getline(inputFile, particleLine)) {
        if (particleLine == "/") {
            // End of particle definition
            break;
        }
        
        std::istringstream lineIss(particleLine);
        std::string paramKey;
        std::string paramValue;
        std::string val1;
        std::string val2;
        
        // NSLog(@"ParamValue: %@", [NSString stringWithUTF8String:particleLine.c_str()]);
        if (lineIss >> paramKey >> paramValue) {
            if (paramKey == "InitialParticleNumber") {
                particle.pNum = std::stoi(paramValue);
            } else if (paramKey == "MaxParticleNumber") {
                particle.pNumMax = std::stoi(paramValue);
            } else if (paramKey == "Charge") {
                particle.q = std::stod(paramValue);
            } else if (paramKey == "Mass") {
                particle.m = std::stod(paramValue);
            } else if (paramKey == "Weight") {
                particle.w = std::stod(paramValue);
            } else if (paramKey == "GenerateType") {
                particle.GenerateType = [NSString stringWithUTF8String:paramValue.c_str()];
            } else if (paramKey == "InitialPosX") {
                particle.initX[0] = std::stod(paramValue);
                if (lineIss >> val1) {
                particle.initX[1] = std::stod(val1);
                }
            } else if (paramKey == "InitialPosY") {
                particle.initY[0] = std::stod(paramValue);
                if (lineIss >> val1) {
                particle.initY[1] = std::stod(val1);
                }
            } else if (paramKey == "InitialVel") {
                particle.initU[0] = std::stod(paramValue);
                if (lineIss >> val1 >> val2) {
                particle.initU[1] = std::stod(val1);
                particle.initU[2] = std::stod(val2);
                }
            } else if (paramKey == "InitialTemp[eV]") {
                particle.initT = std::stod(paramValue)*evtok;
            }
        }
    }    
    // Add the completed particle to our vector
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
    std::string typeLine;
    if (std::getline(inputFile, typeLine)) {
        std::istringstream typeIss(typeLine);
        if (typeIss >> key >> value && key == "BCType") {
            boundary.type = [NSString stringWithUTF8String:value.c_str()];;
        }
    }
    
    // Parse optional value line (if any before the '/' delimiter)
    std::string valueLine;
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
    
    // Add the boundary to our vector
    _particleBoundaries.push_back(boundary);
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
        } else if (key == "InitializeType") {
            _field.InitType = [NSString stringWithUTF8String:value.c_str()];
        } else if (key == "AmplitudeOfE[SI]") {
            _field.ampE[0] = stod(value)*VtoG;
            if (iss >> val1 >> val2) {
                _field.ampE[1] = stod(val1)*VtoG;
                _field.ampE[2] = stod(val2)*VtoG;
            }
        } else if (key == "AmplitudeOfB[SI]") {
            _field.ampB[0] = stod(value)*TtoG;
            if (iss >> val1 >> val2) {
                _field.ampB[1] = stod(val1)*TtoG;
                _field.ampB[2] = stod(val2)*TtoG;
            }
        } else if (key == "WeightingOrder") {
            _field.weightOrder = stoi(value);
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
    std::string typeLine;
    if (std::getline(inputFile, typeLine)) {
        std::istringstream typeIss(typeLine);
        if (typeIss >> key >> value && key == "BCType") {
            boundary.type = [NSString stringWithUTF8String:value.c_str()];
        }
    }
    
    // Parse optional value line (if any before the '/' delimiter)
    std::string valueLine;
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
    
    // Add the boundary to our vector
    _fieldBoundaries.push_back(boundary);
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